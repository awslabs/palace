// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <cmath>
#include <complex>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/qfunctions/coeff/coeff_qf.h"
#include "fem/qfunctions/coeff/pml_qf.h"
#include "linalg/vector.hpp"
#include "models/pml.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/jsonschema.hpp"
#include "utils/units.hpp"

namespace palace
{
using json = nlohmann::json;
using namespace Catch;
using namespace std::complex_literals;

namespace
{

// Cube [-1, 1]³ with boundary attribute 1, optionally with domain attributes 1 for z < -0.5
// and 2 for z > -0.5 (for n a multiple of 4), and nodes of the given order.
mfem::Mesh MakeBox(mfem::Element::Type type, int n, int order, bool layered = false)
{
  auto mesh = mfem::Mesh::MakeCartesian3D(n, n, n, type, 2.0, 2.0, 2.0);
  mesh.Transform(
      [](const mfem::Vector &x, mfem::Vector &y)
      {
        y = x;
        y -= 1.0;
      });
  mesh.SetCurvature(order);
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    mfem::Vector c;
    mesh.GetElementCenter(e, c);
    mesh.SetAttribute(e, (layered && c(2) > -0.5) ? 2 : 1);
  }
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    mesh.SetBdrAttribute(be, 1);
  }
  mesh.SetAttributes();
  return mesh;
}

// Mesh with the vertices at the positions of the nodes.
struct NodalVerticesMesh : public mfem::Mesh
{
  NodalVerticesMesh(mfem::Mesh &&mesh) : mfem::Mesh(std::move(mesh))
  {
    SetVerticesFromNodes(GetNodes());
  }
};

// Unit ball meshed from the cube [-1, 1]³, with the nodes of its boundary on the unit
// sphere (by the smooth map x_i ↦ x_i (1 - x_j² / 2 - x_k² / 2 + x_j² x_k² / 3)^½, which
// maps the faces of the cube onto the sphere).
mfem::Mesh MakeBall(mfem::Element::Type type, int n, int order)
{
  auto mesh = MakeBox(type, n, order);
  mesh.Transform(
      [](const mfem::Vector &x, mfem::Vector &y)
      {
        y.SetSize(3);
        for (int i = 0; i < 3; i++)
        {
          const double a = x((i + 1) % 3) * x((i + 1) % 3),
                       b = x((i + 2) % 3) * x((i + 2) % 3);
          y(i) = x(i) * std::sqrt(1.0 - 0.5 * a - 0.5 * b + a * b / 3.0);
        }
      });
  return NodalVerticesMesh(std::move(mesh));
}

config::PMLData MakePML(json extra = json::object())
{
  json pml = {{"Attributes", {1}}, {"Thickness", 0.5}, {"Layers", 2}};
  pml.update(extra);
  return config::PMLData(pml, true);
}

// Extrude the entry pml of config["Boundaries"]["PML"], whose surface must be one part, and
// replace it with its block.
void ExtrudeEntry(std::unique_ptr<mfem::Mesh> &mesh, config::PMLData &pml,
                  const config::BoundaryData &boundaries = {})
{
  config::DomainData domains;
  domains.pml = {pml};
  mesh::ExtrudePML(mesh, domains, boundaries);
  REQUIRE(domains.pml.size() == 1);
  pml = std::move(domains.pml[0]);
}

// Position of the nodal degree of freedom j of element e (or boundary element) of a mesh
// with H1 nodes.
mfem::Vector NodePosition(const mfem::Mesh &mesh, int e, int j, bool bdr = false)
{
  const auto &nodes = *mesh.GetNodes();
  mfem::Array<int> vdofs;
  bdr ? nodes.FESpace()->GetBdrElementVDofs(e, vdofs)
      : nodes.FESpace()->GetElementVDofs(e, vdofs);
  const int nd = vdofs.Size() / 3;
  mfem::Vector x(3);
  for (int c = 0; c < 3; c++)
  {
    x(c) = nodes(vdofs[j + c * nd]);
  }
  return x;
}

// Smallest determinant of the Jacobian of the extruded elements (from ne) at the vertices
// and quadrature points.
double MinExtrudedJacobian(mfem::Mesh &mesh, int ne, int order)
{
  double det_min = mfem::infinity();
  mfem::IsoparametricTransformation T;
  for (int e = ne; e < mesh.GetNE(); e++)
  {
    mesh.GetElementTransformation(e, &T);
    const auto &ir = mfem::IntRules.Get(T.GetGeometryType(), 2 * order + 2);
    for (int j = 0; j < ir.GetNPoints(); j++)
    {
      T.SetIntPoint(&ir.IntPoint(j));
      det_min = std::min(det_min, T.Jacobian().Det());
    }
  }
  return det_min;
}

}  // namespace

TEST_CASE("ExtrudePML on a sphere", "[pml][geodata][Serial]")
{
  // The layer extruded from a curved sphere is the spherical shell: every node, of every
  // order, is at the radius 1 + ξ for its depth ξ = (k + z) t / N into the layer. (The
  // flat faces of order 1 are too far from the sphere for this layer.)
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  const int order = GENERATE(2, 3);
  auto mesh = std::make_unique<mfem::Mesh>(MakeBall(type, 2, order));
  const int ne = mesh->GetNE(), nbe = mesh->GetNBE();
  auto pml = MakePML();
  ExtrudeEntry(mesh, pml, {});

  // The default core is the centre of the bounding box of the sphere.
  CHECK(pml.extrusion.center[0] == pml.extrusion.center[1]);
  for (int c = 0; c < 3; c++)
  {
    CHECK_THAT(pml.extrusion.center[0][c], Matchers::WithinAbs(0.0, 1.0e-12));
  }
  CHECK(pml.attributes == std::vector<int>{2});
  CHECK(pml.extrusion.base_attributes == std::vector<int>{1});
  CHECK(pml.extrusion.back_attribute == 2);
  REQUIRE(mesh->GetNE() == ne + 2 * nbe);

  double err = 0.0;
  const auto &fespace = *mesh->GetNodes()->FESpace();
  for (int e = ne; e < mesh->GetNE(); e++)
  {
    CHECK(mesh->GetElementGeometry(e) == ((type == mfem::Element::TETRAHEDRON)
                                              ? mfem::Geometry::PRISM
                                              : mfem::Geometry::CUBE));
    CHECK(mesh->GetAttribute(e) == 2);
    const auto &ir = fespace.GetFE(e)->GetNodes();
    for (int j = 0; j < ir.GetNPoints(); j++)
    {
      const double xi = ((e - ne) % 2 + ir.IntPoint(j).z) * 0.25;
      err = std::max(err, std::abs(NodePosition(*mesh, e, j).Norml2() - (1.0 + xi)));
    }
  }
  CHECK(err < 1.0e-12);

  // The sphere is an interior boundary, and the back is the sphere of radius 1.5.
  int n_back = 0;
  for (int be = 0; be < mesh->GetNBE(); be++)
  {
    int e1, e2;
    mesh->GetFaceElements(mesh->GetBdrElementFaceIndex(be), &e1, &e2);
    if (mesh->GetBdrAttribute(be) == 1)
    {
      CHECK(e2 >= 0);
      continue;
    }
    REQUIRE(mesh->GetBdrAttribute(be) == 2);
    CHECK(e2 < 0);
    n_back++;
    const int nd = fespace.GetBE(be)->GetDof();
    for (int j = 0; j < nd; j++)
    {
      CHECK_THAT(NodePosition(*mesh, be, j, true).Norml2(),
                 Matchers::WithinAbs(1.5, 1.0e-12));
    }
  }
  CHECK(n_back == nbe);
}

TEST_CASE("ExtrudePML on a box", "[pml][geodata][Serial]")
{
  // Columns of elements at the edges and corners of a box are valid and conforming: the
  // extruded elements have positive Jacobians, and the only exterior faces of the mesh are
  // those of the back of the layer, whose nodes are on the box grown by the thickness.
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  const int order = GENERATE(1, 2);
  auto mesh = std::make_unique<mfem::Mesh>(MakeBox(type, 2, order));
  const int ne = mesh->GetNE(), nbe = mesh->GetNBE();
  auto pml = MakePML();
  ExtrudeEntry(mesh, pml, {});
  REQUIRE(mesh->GetNE() == ne + 2 * nbe);
  CHECK(MinExtrudedJacobian(*mesh, ne, order) > 0.0);
  int n_exterior = 0;
  for (int f = 0; f < mesh->GetNumFaces(); f++)
  {
    int e1, e2;
    mesh->GetFaceElements(f, &e1, &e2);
    n_exterior += (e2 < 0);
  }
  CHECK(n_exterior == nbe);
  for (int be = 0; be < mesh->GetNBE(); be++)
  {
    if (mesh->GetBdrAttribute(be) != pml.extrusion.back_attribute)
    {
      continue;
    }
    for (int j = 0; j < mesh->GetNodes()->FESpace()->GetBE(be)->GetDof(); j++)
    {
      CHECK_THAT(NodePosition(*mesh, be, j, true).Normlinf(),
                 Matchers::WithinAbs(1.5, 1.0e-12));
    }
  }
}

TEST_CASE("ExtrudePML numbers the new attributes above the configuration",
          "[pml][geodata][Serial]")
{
  // The box has the domain attribute 1 and the boundary attribute 1. A material for the
  // domain attribute 2 and a boundary condition for 5 and a postprocessing surface for 6
  // (or the other way round), none of which is in the mesh, must not apply to the layer or
  // its back.
  const bool swap = GENERATE(false, true);
  auto mesh = std::make_unique<mfem::Mesh>(MakeBox(mfem::Element::HEXAHEDRON, 2, 1));
  config::DomainData domains;
  domains.attributes = {1, 2};
  domains.pml = {MakePML()};
  config::BoundaryData boundaries;
  boundaries.attributes = {swap ? 6 : 5};
  boundaries.postpro.attributes = {swap ? 5 : 6};
  mesh::ExtrudePML(mesh, domains, boundaries);
  REQUIRE(domains.pml.size() == 1);
  CHECK(domains.pml[0].attributes == std::vector<int>{3});
  CHECK(domains.pml[0].extrusion.back_attribute == 7);
}

TEST_CASE("ExtrudePML with a segment core", "[pml][geodata][Serial]")
{
  // With a segment core through the interface z = -0.5 between two materials of the box,
  // the nodes beside the segment move normal to it, and the interface stays planar through
  // the layer.
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  auto mesh = std::make_unique<mfem::Mesh>(MakeBox(type, 4, 2, true));
  const int ne = mesh->GetNE();
  auto pml = MakePML({{"Center", {{0.0, 0.0, -0.75}, {0.0, 0.0, 0.75}}}});
  ExtrudeEntry(mesh, pml, {});
  CHECK(pml.attributes == std::vector<int>{3, 4});
  CHECK(pml.extrusion.base_attributes == std::vector<int>{1, 2});
  CHECK(MinExtrudedJacobian(*mesh, ne, 2) > 0.0);
  double z_below = -mfem::infinity(), z_above = mfem::infinity();
  for (int e = ne; e < mesh->GetNE(); e++)
  {
    for (int j = 0; j < mesh->GetNodes()->FESpace()->GetFE(e)->GetDof(); j++)
    {
      const double z = NodePosition(*mesh, e, j)(2);
      (mesh->GetAttribute(e) == 3) ? (z_below = std::max(z_below, z))
                                   : (z_above = std::min(z_above, z));
    }
  }
  CHECK_THAT(z_below, Matchers::WithinAbs(-0.5, 1.0e-12));
  CHECK_THAT(z_above, Matchers::WithinAbs(-0.5, 1.0e-12));

  // With the default point core, the interface would tilt into a cone, and the error gives
  // the core grown to span it, which keeps it flat.
  mesh = std::make_unique<mfem::Mesh>(MakeBox(type, 4, 2, true));
  pml = MakePML();
  CHECK_THROWS_WITH(
      ExtrudeEntry(mesh, pml, {}),
      Matchers::ContainsSubstring("The interface between the materials of the domain "
                                  "attributes 1 and 2 crosses") &&
          Matchers::ContainsSubstring("\"Center\": [[0, 0, -0.5], [0, 0, 0]]"));
  mesh = std::make_unique<mfem::Mesh>(MakeBox(type, 4, 2, true));
  pml = MakePML({{"Center", {{0.0, 0.0, -0.5}, {0.0, 0.0, 0.0}}}});
  CHECK_NOTHROW(ExtrudeEntry(mesh, pml, {}));
}

TEST_CASE("ExtrudePML keeps a boundary across a box flat", "[pml][geodata][Serial]")
{
  // An interior PEC boundary (attribute 2) of the box, of the faces with all of their
  // vertices on it, through the PML (attribute 1): on the plane z = z₀, a strip |y| < 0.5
  // for x > 0 through the face x = 1, or a sheet through the whole box; on the plane x = y,
  // which runs into the vertical edges of the box, a sheet through the whole box, or a fin
  // x > 0.5, |z| < 0.5, which touches the box only along its edge.
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  auto MakeMesh = [&](auto &&on)
  {
    auto box = MakeBox(type, 4, 1);
    std::vector<int> sheet;
    mfem::Array<int> verts;
    for (int f = 0; f < box.GetNumFaces(); f++)
    {
      int e1, e2;
      box.GetFaceElements(f, &e1, &e2);
      box.GetFaceVertices(f, verts);
      bool all = (e2 >= 0);
      for (const int v : verts)
      {
        all = all && on(box.GetVertex(v));
      }
      if (all)
      {
        sheet.push_back(f);
      }
    }
    auto mesh = std::make_unique<mfem::Mesh>(
        3, box.GetNV(), box.GetNE(), box.GetNBE() + static_cast<int>(sheet.size()), 3);
    for (int v = 0; v < box.GetNV(); v++)
    {
      mesh->AddVertex(box.GetVertex(v));
    }
    for (int e = 0; e < box.GetNE(); e++)
    {
      mesh->AddElement(box.GetElement(e)->Duplicate(mesh.get()));
    }
    for (int be = 0; be < box.GetNBE(); be++)
    {
      mesh->AddBdrElement(box.GetBdrElement(be)->Duplicate(mesh.get()));
    }
    for (const int f : sheet)
    {
      auto *el = box.GetFace(f)->Duplicate(mesh.get());
      el->SetAttribute(2);
      mesh->AddBdrElement(el);
    }
    mesh->FinalizeTopology(false);
    mesh->SetCurvature(1);
    return mesh;
  };
  auto Strip = [](const double *x)
  { return std::abs(x[2]) < 1.0e-12 && x[0] > -1.0e-12 && std::abs(x[1]) < 0.5 + 1.0e-12; };
  auto Sheet = [](double z0)
  { return [z0](const double *x) { return std::abs(x[2] - z0) < 1.0e-12; }; };
  config::BoundaryData pec;
  pec.attributes = {2};

  // The default core, the centre of the box, would widen the strip in the layer: the error
  // gives the core grown to span it.
  auto mesh = MakeMesh(Strip);
  auto pml = MakePML();
  CHECK_THROWS_WITH(
      ExtrudeEntry(mesh, pml, pec),
      Matchers::ContainsSubstring("Boundary attribute 2 crosses") &&
          Matchers::ContainsSubstring(" at (1, ") &&
          Matchers::ContainsSubstring("\"Center\": [[0, -0.5, 0], [0, 0.5, 0]]"));

  // A boundary without a boundary condition does not need to stay flat.
  mesh = MakeMesh(Strip);
  pml = MakePML();
  CHECK_NOTHROW(ExtrudeEntry(mesh, pml, {}));

  // With that core, the strip continues straight and flat through the layer, with one
  // quadrilateral per layer on each of its two edges on the face.
  auto CountStrip = [](const mfem::Mesh &mesh)
  {
    int count = 0;
    for (int be = 0; be < mesh.GetNBE(); be++)
    {
      count += (mesh.GetBdrAttribute(be) == 2);
    }
    return count;
  };
  // The largest |x_c| and |x_0 - x_1| of the nodes of the boundary.
  auto Extent = [](const mfem::Mesh &mesh)
  {
    std::array<double, 4> ext{};
    for (int be = 0; be < mesh.GetNBE(); be++)
    {
      for (int j = 0; mesh.GetBdrAttribute(be) == 2 &&
                      j < mesh.GetNodes()->FESpace()->GetBE(be)->GetDof();
           j++)
      {
        const auto x = NodePosition(mesh, be, j, true);
        for (int c = 0; c < 3; c++)
        {
          ext[c] = std::max(ext[c], std::abs(x(c)));
        }
        ext[3] = std::max(ext[3], std::abs(x(0) - x(1)));
      }
    }
    return ext;
  };
  mesh = MakeMesh(Strip);
  const int n_strip = CountStrip(*mesh);
  pml = MakePML({{"Center", {{0.0, -0.5, 0.0}, {0.0, 0.5, 0.0}}}});
  ExtrudeEntry(mesh, pml, pec);
  CHECK(CountStrip(*mesh) == n_strip + 2 * 2);
  CHECK_THAT(Extent(*mesh)[1], Matchers::WithinAbs(0.5, 1.0e-12));
  CHECK(Extent(*mesh)[2] < 1.0e-12);

  // The sheet through the whole box stays flat with the default core, which is on it, and
  // not on the plane z = 0.5, off the core.
  mesh = MakeMesh(Sheet(0.0));
  pml = MakePML();
  CHECK_NOTHROW(ExtrudeEntry(mesh, pml, pec));
  mesh = MakeMesh(Sheet(0.5));
  pml = MakePML();
  CHECK_THROWS_WITH(
      ExtrudeEntry(mesh, pml, pec),
      Matchers::ContainsSubstring("Boundary attribute 2 crosses") &&
          Matchers::ContainsSubstring("\"Center\": [[0, 0, 0], [0, 0, 0.5]]"));
  mesh = MakeMesh(Sheet(0.5));
  pml = MakePML({{"Center", {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.5}}}});
  CHECK_NOTHROW(ExtrudeEntry(mesh, pml, pec));

  // A boundary along an edge of the box continues there. On the plane x = y, which has
  // faces on the tetrahedral mesh only, the sheet stays in its plane with the default core,
  // out to the back of the layer, with its 16 edges on the box (4 on the top and the bottom
  // and on each of the two edges of the box) extended. The fin, with its two edges along
  // the edge of the box, would grow in the layer: the error gives the core which spans it.
  if (type != mfem::Element::TETRAHEDRON)
  {
    return;
  }
  auto Diagonal = [](const double *x) { return std::abs(x[0] - x[1]) < 1.0e-12; };
  mesh = MakeMesh(Diagonal);
  const int n_sheet = CountStrip(*mesh);
  REQUIRE(n_sheet == 2 * 4 * 4);
  pml = MakePML();
  CHECK_NOTHROW(ExtrudeEntry(mesh, pml, pec));
  CHECK(CountStrip(*mesh) == n_sheet + 16 * 2);
  CHECK_THAT(Extent(*mesh)[0], Matchers::WithinAbs(1.5, 1.0e-12));
  CHECK(Extent(*mesh)[3] < 1.0e-12);
  auto Fin = [](const double *x)
  {
    return std::abs(x[0] - x[1]) < 1.0e-12 && x[0] > 0.5 - 1.0e-12 &&
           std::abs(x[2]) < 0.5 + 1.0e-12;
  };
  mesh = MakeMesh(Fin);
  const int n_fin = CountStrip(*mesh);
  REQUIRE(n_fin == 2 * 2);
  pml = MakePML();
  CHECK_THROWS_WITH(
      ExtrudeEntry(mesh, pml, pec),
      Matchers::ContainsSubstring("Boundary attribute 2 crosses") &&
          Matchers::ContainsSubstring(" at (1, 1, ") &&
          Matchers::ContainsSubstring("\"Center\": [[0, 0, -0.5], [0, 0, 0.5]]"));
  mesh = MakeMesh(Fin);
  pml = MakePML({{"Center", {{0.0, 0.0, -0.5}, {0.0, 0.0, 0.5}}}});
  ExtrudeEntry(mesh, pml, pec);
  CHECK(CountStrip(*mesh) == n_fin + 2 * 2);
  CHECK_THAT(Extent(*mesh)[2], Matchers::WithinAbs(0.5, 1.0e-12));
  CHECK(Extent(*mesh)[3] < 1.0e-12);
}

TEST_CASE("ExtrudePML rejects invalid surfaces", "[pml][geodata][Serial]")
{
  using Catch::Matchers::ContainsSubstring;
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  auto Extrude =
      [&](json extra, const config::BoundaryData &boundaries = {}, bool interior = false)
  {
    auto mesh = std::make_unique<mfem::Mesh>(MakeBox(type, 2, 1, interior));
    if (interior)
    {
      // Boundary elements with attribute 1 on the interior faces z = 0 as well.
      std::vector<int> interior_faces;
      for (int f = 0; f < mesh->GetNumFaces(); f++)
      {
        int e1, e2;
        mesh->GetFaceElements(f, &e1, &e2);
        if (e2 >= 0 && mesh->GetAttribute(e1) != mesh->GetAttribute(e2))
        {
          interior_faces.push_back(f);
        }
      }
      auto copy = std::make_unique<mfem::Mesh>(
          3, mesh->GetNV(), mesh->GetNE(),
          mesh->GetNBE() + static_cast<int>(interior_faces.size()), 3);
      for (int v = 0; v < mesh->GetNV(); v++)
      {
        copy->AddVertex(mesh->GetVertex(v));
      }
      for (int e = 0; e < mesh->GetNE(); e++)
      {
        copy->AddElement(mesh->GetElement(e)->Duplicate(copy.get()));
      }
      for (int be = 0; be < mesh->GetNBE(); be++)
      {
        copy->AddBdrElement(mesh->GetBdrElement(be)->Duplicate(copy.get()));
      }
      for (const int f : interior_faces)
      {
        auto *el = mesh->GetFace(f)->Duplicate(copy.get());
        el->SetAttribute(1);
        copy->AddBdrElement(el);
      }
      copy->FinalizeTopology(false);
      copy->SetCurvature(1);
      mesh = std::move(copy);
    }
    auto pml = MakePML(extra);
    ExtrudeEntry(mesh, pml, boundaries);
  };

  // A core outside of the box or on its surface, which is not inside of its side x = 1.
  CHECK_THROWS_WITH(Extrude({{"Center", {3.0, 0.0, 0.0}}}),
                    ContainsSubstring("is not inside of its side"));
  CHECK_THROWS_WITH(Extrude({{"Center", {1.0, 0.0, 0.0}}}),
                    ContainsSubstring("is not inside of its side"));
  // An interior boundary.
  CHECK_THROWS_WITH(Extrude(json::object(), {}, true),
                    ContainsSubstring("is not on the exterior boundary"));

  // A flat face of the box, along its normal, or along the rays from a core behind it.
  auto ExtrudeFace = [&](json extra)
  {
    auto mesh = std::make_unique<mfem::Mesh>(MakeBox(type, 2, 1));
    for (int be = 0; be < mesh->GetNBE(); be++)
    {
      mfem::Vector c(3);
      mfem::IsoparametricTransformation T;
      mesh->GetBdrElementTransformation(be, &T);
      T.Transform(mfem::Geometries.GetCenter(T.GetGeometryType()), c);
      mesh->SetBdrAttribute(be, (c(2) > 1.0 - 1.0e-12) ? 1 : 2);
    }
    mesh->SetAttributes();
    auto pml = MakePML(extra);
    ExtrudeEntry(mesh, pml, {});
    return pml.extrusion.depth;
  };
  CHECK(ExtrudeFace(json::object()) == config::PMLData::Extrusion::Depth::NORMAL);
  // "Center" only shapes the parts which are not flat, and the radial extrusion rejects a
  // flat face.
  CHECK_THROWS_WITH(ExtrudeFace({{"Center", {0.0, 0.0, -10.0}}}),
                    ContainsSubstring("\"Center\" is not used by the PML boundary with "
                                      "attributes 1, whose parts are flat"));
  CHECK_THROWS_WITH(ExtrudeFace({{"Extrusion", "Radial"}}),
                    ContainsSubstring("which is flat! A flat face is extruded along its "
                                      "normal"));
  // A sphere of flat faces, whose distance from its centre varies by
  // more than the tolerance for a thin layer.
  {
    auto mesh = std::make_unique<mfem::Mesh>(MakeBall(type, 2, 1));
    auto pml = MakePML({{"Thickness", 0.01}});
    CHECK_THROWS_WITH(ExtrudeEntry(mesh, pml, {}),
                      ContainsSubstring("The extruded PML does not fit") &&
                          ContainsSubstring("a sphere around its centre"));
    mesh = std::make_unique<mfem::Mesh>(MakeBall(type, 2, 1));
    pml = MakePML({{"Thickness", 10.0}});
    CHECK_NOTHROW(ExtrudeEntry(mesh, pml, {}));
    CHECK(pml.extrusion.depth == config::PMLData::Extrusion::Depth::DISTANCE);
  }
  // A surface which fits none of the shapes: an L-shaped box, whose faces are not all on
  // its bounding box, and a cube with faces normal to no axis. Along their normals, the
  // faces of the L-shaped box would share their edges.
  for (const bool normal : {false, true})
  {
    auto cubes = mfem::Mesh::MakeCartesian3D(2, 2, 1, type, 2.0, 2.0, 1.0);
    std::vector<int> keep;
    for (int e = 0; e < cubes.GetNE(); e++)
    {
      mfem::Vector c;
      cubes.GetElementCenter(e, c);
      if (c(0) < 1.0 || c(1) < 1.0)
      {
        keep.push_back(e);
      }
    }
    auto mesh =
        std::make_unique<mfem::Mesh>(3, cubes.GetNV(), static_cast<int>(keep.size()), 0, 3);
    for (int v = 0; v < cubes.GetNV(); v++)
    {
      mesh->AddVertex(cubes.GetVertex(v));
    }
    for (const int e : keep)
    {
      mesh->AddElement(cubes.GetElement(e)->Duplicate(mesh.get()));
    }
    mesh->FinalizeTopology(true);
    mesh->RemoveUnusedVertices();
    for (int be = 0; be < mesh->GetNBE(); be++)
    {
      mesh->SetBdrAttribute(be, 1);
    }
    mesh->SetAttributes();
    mesh->SetCurvature(1);
    auto pml = normal ? MakePML({{"Extrusion", "Normal"}}) : MakePML();
    CHECK_THROWS_WITH(ExtrudeEntry(mesh, pml, {}),
                      normal
                          ? ContainsSubstring("which share the edge at") &&
                                ContainsSubstring("is for separate flat faces")
                          : ContainsSubstring("is not on the bounding box of the part") &&
                                !ContainsSubstring("\"Normal\""));
  }
  {
    auto mesh = std::make_unique<mfem::Mesh>(MakeBox(type, 2, 1));
    mesh->Transform(
        [](const mfem::Vector &x, mfem::Vector &y)
        {
          y.SetSize(3);
          const double c = std::cos(0.3), s = std::sin(0.3);
          y(0) = c * x(0) - s * x(1);
          y(1) = s * x(0) + c * x(1);
          y(2) = x(2);
        });
    auto pml = MakePML();
    CHECK_THROWS_WITH(ExtrudeEntry(mesh, pml, {}),
                      ContainsSubstring("or a box with its faces normal to the coordinate "
                                        "axes") &&
                          !ContainsSubstring("\"Normal\""));
  }

  // A half sphere on a ground plane (attribute 2) ends flush with it for its centre on the
  // plane, and not for a centre above it.
  auto ExtrudeHalfBall = [&](json extra)
  {
    auto half = mfem::Mesh::MakeCartesian3D(2, 2, 1, type, 2.0, 2.0, 1.0);
    half.SetCurvature(2);
    half.Transform(
        [](const mfem::Vector &x, mfem::Vector &y)
        {
          const double z[3] = {x(0) - 1.0, x(1) - 1.0, x(2)};
          y.SetSize(3);
          for (int i = 0; i < 3; i++)
          {
            const double a = z[(i + 1) % 3] * z[(i + 1) % 3],
                         b = z[(i + 2) % 3] * z[(i + 2) % 3];
            y(i) = z[i] * std::sqrt(1.0 - 0.5 * a - 0.5 * b + a * b / 3.0);
          }
        });
    for (int be = 0; be < half.GetNBE(); be++)
    {
      half.SetBdrAttribute(be, (half.GetBdrAttribute(be) == 1) ? 2 : 1);
    }
    half.SetAttributes();
    auto mesh = std::make_unique<mfem::Mesh>(NodalVerticesMesh(std::move(half)));
    json entry = {{"Thickness", 2.0}};
    entry.update(extra);
    auto pml = MakePML(entry);
    ExtrudeEntry(mesh, pml, {});
    return pml.extrusion.depth;
  };
  CHECK(ExtrudeHalfBall({{"Center", {0.0, 0.0, 0.0}}}) ==
        config::PMLData::Extrusion::Depth::DISTANCE);
  CHECK_THROWS_WITH(ExtrudeHalfBall({{"Center", {0.0, 0.0, 0.05}}}),
                    ContainsSubstring("does not end flush with the boundary attribute 2"));

  // No boundary elements, and a back attribute of the mesh.
  CHECK_THROWS_WITH(Extrude({{"Attributes", {7}}}),
                    ContainsSubstring("No boundary elements"));
  CHECK_THROWS_WITH(Extrude({{"BackAttribute", 1}}),
                    ContainsSubstring("is already a boundary attribute"));

  // Ports, periodic, and Floquet boundaries which touch the surface: the surface is the box
  // except its face x = -1 (attribute 2).
  auto ExtrudeWithPort = [&](const config::BoundaryData &boundaries)
  {
    auto mesh = std::make_unique<mfem::Mesh>(MakeBox(type, 2, 1));
    for (int be = 0; be < mesh->GetNBE(); be++)
    {
      mfem::Vector c(3);
      mfem::IsoparametricTransformation T;
      mesh->GetBdrElementTransformation(be, &T);
      T.Transform(mfem::Geometries.GetCenter(T.GetGeometryType()), c);
      mesh->SetBdrAttribute(be, (c(0) < -1.0 + 1.0e-12) ? 2 : 1);
    }
    mesh->SetAttributes();
    auto pml = MakePML();
    ExtrudeEntry(mesh, pml, boundaries);
  };
  config::BoundaryData boundaries;
  CHECK_NOTHROW(ExtrudeWithPort(boundaries));
  boundaries.lumpedport[1].elements.emplace_back().attributes = {2};
  CHECK_THROWS_WITH(ExtrudeWithPort(boundaries),
                    ContainsSubstring("touches the PML boundary"));
  boundaries.lumpedport.clear();
  boundaries.waveport[1].attributes = {2};
  CHECK_THROWS_WITH(ExtrudeWithPort(boundaries),
                    ContainsSubstring("touches the PML boundary"));
  boundaries.waveport.clear();
  boundaries.periodic.boundary_pairs.emplace_back().donor_attributes = {2};
  CHECK_THROWS_WITH(ExtrudeWithPort(boundaries),
                    ContainsSubstring("touches the PML boundary"));
}

TEST_CASE("ExtrudePML extends other boundaries through the layer", "[pml][geodata][Serial]")
{
  // Half of the box above a ground plane z = -1 (attribute 2), which the extrusion extends
  // into the quadrilaterals below the layer of its edges with the box.
  auto mesh = std::make_unique<mfem::Mesh>(MakeBox(mfem::Element::TETRAHEDRON, 2, 2));
  for (int be = 0; be < mesh->GetNBE(); be++)
  {
    mfem::Vector c(3);
    mfem::IsoparametricTransformation T;
    mesh->GetBdrElementTransformation(be, &T);
    T.Transform(mfem::Geometries.GetCenter(T.GetGeometryType()), c);
    mesh->SetBdrAttribute(be, (c(2) < -1.0 + 1.0e-12) ? 2 : 1);
  }
  mesh->SetAttributes();
  const int nbe = mesh->GetNBE();
  auto pml = MakePML({{"Center", {{0.0, 0.0, -1.0}, {0.0, 0.0, 0.5}}}});
  ExtrudeEntry(mesh, pml, {});
  // The ground plane has 8 boundary edges with the box, each extended into 2 quadrilaterals
  // in the plane z = -1, and the mesh has no exterior face without a boundary element.
  int n_ground = 0;
  for (int be = nbe; be < mesh->GetNBE(); be++)
  {
    if (mesh->GetBdrAttribute(be) == 2)
    {
      n_ground++;
      CHECK(mesh->GetBdrElementGeometry(be) == mfem::Geometry::SQUARE);
      for (int j = 0; j < mesh->GetNodes()->FESpace()->GetBE(be)->GetDof(); j++)
      {
        CHECK_THAT(NodePosition(*mesh, be, j, true)(2), Matchers::WithinAbs(-1.0, 1.0e-12));
      }
    }
  }
  CHECK(n_ground == 16);
  std::vector<bool> has_bdr(mesh->GetNumFaces(), false);
  for (int be = 0; be < mesh->GetNBE(); be++)
  {
    has_bdr[mesh->GetBdrElementFaceIndex(be)] = true;
  }
  for (int f = 0; f < mesh->GetNumFaces(); f++)
  {
    int e1, e2;
    mesh->GetFaceElements(f, &e1, &e2);
    CHECK((e2 >= 0 || has_bdr[f]));
  }
}

namespace
{

// Box [-1, 1]³ (see MakeBox) with the boundary attribute 2 on its face z = 1.
mfem::Mesh MakeBoxTop(mfem::Element::Type type, int n, int order)
{
  auto mesh = MakeBox(type, n, order);
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    mfem::Vector c(3);
    mfem::IsoparametricTransformation T;
    mesh.GetBdrElementTransformation(be, &T);
    T.Transform(mfem::Geometries.GetCenter(T.GetGeometryType()), c);
    mesh.SetBdrAttribute(be, (c(2) > 1.0 - 1.0e-12) ? 2 : 1);
  }
  mesh.SetAttributes();
  return mesh;
}

// Cylinder of radius 1 around the z-axis, z ∈ [-1, 1], meshed from the box [-1, 1]³, with
// the nodes of its wall (boundary attribute 1) on the cylinder (by the smooth map (x, y) ↦
// (x (1 - y² / 2)^½, y (1 - x² / 2)^½) of the square onto the disk), and its caps of
// boundary attribute 2.
mfem::Mesh MakeCylinder(mfem::Element::Type type, int n, int order)
{
  auto mesh = MakeBox(type, n, order);
  mesh.Transform(
      [](const mfem::Vector &x, mfem::Vector &y)
      {
        y = x;
        y(0) = x(0) * std::sqrt(1.0 - 0.5 * x(1) * x(1));
        y(1) = x(1) * std::sqrt(1.0 - 0.5 * x(0) * x(0));
      });
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    mfem::Vector c(3);
    mfem::IsoparametricTransformation T;
    mesh.GetBdrElementTransformation(be, &T);
    T.Transform(mfem::Geometries.GetCenter(T.GetGeometryType()), c);
    mesh.SetBdrAttribute(be, (std::abs(c(2)) > 1.0 - 1.0e-12) ? 2 : 1);
  }
  mesh.SetAttributes();
  return NodalVerticesMesh(std::move(mesh));
}

// Stretch of the PML block of an extruded layer.
pml::Stretch ExtrudedStretch(const config::PMLData &pml)
{
  pml::Stretch p;
  p.extruded = true;
  p.geometry.thickness[0] = pml.extrusion.thickness;
  p.depth = pml.extrusion.depth;
  p.center = pml.extrusion.center;
  p.distance = pml.extrusion.distance;
  p.plane = pml.extrusion.plane;
  return p;
}

// Largest difference over the nodes of the extruded elements of a block of the closed form
// depth and of the given depth of the node positions from the depth (k + z) t / N at which
// the extrusion placed the node, for its layer k and reference coordinate z, or only of the
// closed form from the given depth (on a refined mesh).
template <typename Depth>
double DepthError(const mfem::Mesh &mesh, const config::PMLData &pml, Depth &&depth,
                  bool layers = true)
{
  const auto p = ExtrudedStretch(pml);
  const int N = pml.extrusion.layers;
  double err = 0.0;
  for (int e = 0, i = 0; e < mesh.GetNE(); e++)
  {
    if (!std::ranges::binary_search(pml.attributes, mesh.GetAttribute(e)))
    {
      continue;
    }
    // The extruded elements of a block are in the order of their faces, with the N layers
    // of each face in order.
    const auto &ir = mesh.GetNodes()->FESpace()->GetFE(e)->GetNodes();
    for (int j = 0; j < ir.GetNPoints(); j++)
    {
      const auto x = NodePosition(mesh, e, j);
      const double xi = pml::EvaluateDepth(p, {x(0), x(1), x(2)}).xi, xi_x = depth(x),
                   xi_k = (i % N + ir.IntPoint(j).z) * pml.extrusion.thickness / N;
      err = std::max(err, std::abs(xi - xi_x));
      if (layers)
      {
        err = std::max({err, std::abs(xi - xi_k), std::abs(xi_x - xi_k)});
      }
    }
    i++;
  }
  return err;
}

// Largest difference over the nodes of the extruded elements of a box block of their depth
// beyond the box, the largest distance max_i η_i beyond its PML sides, from the depth (k +
// z) t / N at which the extrusion placed the node: the layer is the box grown by the
// thickness on its PML sides.
double BoxDepthError(const mfem::Mesh &mesh, const config::PMLData &pml)
{
  const auto &ex = pml.extrusion;
  REQUIRE(ex.depth == config::PMLData::Extrusion::Depth::BOX);
  double err = 0.0;
  for (int e = 0, i = 0; e < mesh.GetNE(); e++)
  {
    if (!std::ranges::binary_search(pml.attributes, mesh.GetAttribute(e)))
    {
      continue;
    }
    const auto &ir = mesh.GetNodes()->FESpace()->GetFE(e)->GetNodes();
    for (int j = 0; j < ir.GetNPoints(); j++)
    {
      const auto x = NodePosition(mesh, e, j);
      double eta = 0.0;
      for (int f = 0; f < 6; f++)
      {
        const int a = f / 2;
        eta = ex.box_sides[f]
                  ? std::max(eta, (f % 2) ? x(a) - ex.box[1][a] : ex.box[0][a] - x(a))
                  : eta;
      }
      err = std::max(err, std::abs(eta - (i % ex.layers + ir.IntPoint(j).z) * ex.thickness /
                                             ex.layers));
    }
    i++;
  }
  return err;
}

}  // namespace

namespace
{

// Box [-1, 1]³ (see MakeBox) scaled by s along each axis and shifted by c.
mfem::Mesh MakeScaledBox(mfem::Element::Type type, int n, int order,
                         const std::array<double, 3> &s, const std::array<double, 3> &c,
                         bool layered = false)
{
  auto mesh = MakeBox(type, n, order, layered);
  mesh.Transform(
      [&](const mfem::Vector &x, mfem::Vector &y)
      {
        y.SetSize(3);
        for (int i = 0; i < 3; i++)
        {
          y(i) = c[i] + s[i] * x(i);
        }
      });
  return NodalVerticesMesh(std::move(mesh));
}

// Largest distance of the nodes of the boundary elements with the attribute from the
// surface of the box with the corners lo and hi.
double DistanceFromBox(const mfem::Mesh &mesh, int attr, const std::array<double, 3> &lo,
                       const std::array<double, 3> &hi)
{
  double err = 0.0;
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    for (int j = 0; mesh.GetBdrAttribute(be) == attr &&
                    j < mesh.GetNodes()->FESpace()->GetBE(be)->GetDof();
         j++)
    {
      const auto x = NodePosition(mesh, be, j, true);
      double dist = -mfem::infinity();
      for (int a = 0; a < 3; a++)
      {
        dist = std::max({dist, x(a) - hi[a], lo[a] - x(a)});
      }
      err = std::max(err, std::abs(dist));
    }
  }
  return err;
}

}  // namespace

TEST_CASE("ExtrudePML of a box", "[pml][geodata][Serial]")
{
  // A box which is not a cube is extruded per component from its centre into the box grown
  // by the thickness on every side: the layer has the thickness normal to every face.
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  const int order = GENERATE(1, 2);
  {
    auto mesh = std::make_unique<mfem::Mesh>(
        MakeScaledBox(type, 2, order, {1.0, 0.5, 0.25}, {0.3, 0.0, -0.2}));
    const int ne = mesh->GetNE();
    auto pml = MakePML();
    ExtrudeEntry(mesh, pml, {});
    REQUIRE(pml.extrusion.depth == config::PMLData::Extrusion::Depth::BOX);
    CHECK(BoxDepthError(*mesh, pml) < 1.0e-12);
    CHECK(MinExtrudedJacobian(*mesh, ne, order) > 0.0);
    CHECK(DistanceFromBox(*mesh, pml.extrusion.back_attribute, {-1.2, -1.0, -0.95},
                          {1.8, 1.0, 0.55}) < 1.0e-12);
  }

  // The box of a patch antenna, [-1, 1] x [0, 2] x [0, 1] with a substrate below z = 0.25:
  // PML on the top, the two sides and the back (attribute 7), and open on the port face y =
  // 0 (attribute 8) and the ground plane z = 0 (attribute 9). With the core on both of
  // them, as a segment through the substrate, the layer ends flush with them, their
  // boundaries continue across its ends, and the substrate stays flat.
  auto MakePatch = [&]()
  {
    auto mesh = std::make_unique<mfem::Mesh>(
        MakeScaledBox(type, 4, order, {1.0, 1.0, 0.5}, {0.0, 1.0, 0.5}, true));
    for (int be = 0; be < mesh->GetNBE(); be++)
    {
      mfem::Vector c(3);
      mfem::IsoparametricTransformation T;
      mesh->GetBdrElementTransformation(be, &T);
      T.Transform(mfem::Geometries.GetCenter(T.GetGeometryType()), c);
      mesh->SetBdrAttribute(be, (c(1) < 1.0e-12) ? 8 : (c(2) < 1.0e-12) ? 9 : 7);
    }
    mesh->SetAttributes();
    return mesh;
  };
  auto mesh = MakePatch();
  const int ne = mesh->GetNE(), nbe = mesh->GetNBE();
  auto pml =
      MakePML({{"Attributes", {7}}, {"Center", {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.25}}}});
  ExtrudeEntry(mesh, pml, {});
  REQUIRE(pml.extrusion.depth == config::PMLData::Extrusion::Depth::BOX);
  CHECK(pml.extrusion.box_sides ==
        std::array<bool, 6>{true, true, false, true, false, true});
  CHECK(BoxDepthError(*mesh, pml) < 1.0e-12);
  CHECK(MinExtrudedJacobian(*mesh, ne, order) > 0.0);
  CHECK(DistanceFromBox(*mesh, pml.extrusion.back_attribute, {-1.5, 0.0, 0.0},
                        {1.5, 2.5, 1.5}) < 1.0e-12);
  double y_min = mfem::infinity(), z_min = mfem::infinity(), z_below = -mfem::infinity(),
         z_above = mfem::infinity();
  for (int e = ne; e < mesh->GetNE(); e++)
  {
    for (int j = 0; j < mesh->GetNodes()->FESpace()->GetFE(e)->GetDof(); j++)
    {
      const auto x = NodePosition(*mesh, e, j);
      y_min = std::min(y_min, x(1));
      z_min = std::min(z_min, x(2));
      (mesh->GetAttribute(e) == pml.attributes[0]) ? (z_below = std::max(z_below, x(2)))
                                                   : (z_above = std::min(z_above, x(2)));
    }
  }
  CHECK_THAT(y_min, Matchers::WithinAbs(0.0, 1.0e-12));
  CHECK_THAT(z_min, Matchers::WithinAbs(0.0, 1.0e-12));
  CHECK_THAT(z_below, Matchers::WithinAbs(0.25, 1.0e-12));
  CHECK_THAT(z_above, Matchers::WithinAbs(0.25, 1.0e-12));
  for (const int attr : {8, 9})
  {
    int n_new = 0;
    for (int be = nbe; be < mesh->GetNBE(); be++)
    {
      for (int j = 0; mesh->GetBdrAttribute(be) == attr &&
                      j < mesh->GetNodes()->FESpace()->GetBE(be)->GetDof();
           j++)
      {
        n_new += (j == 0);
        CHECK_THAT(NodePosition(*mesh, be, j, true)((attr == 8) ? 1 : 2),
                   Matchers::WithinAbs(0.0, 1.0e-12));
      }
    }
    CHECK(n_new > 0);
  }

  // The default core is the centre of the box moved onto the open sides, which keeps the
  // substrate flat too, and a core off the open sides is an error.
  mesh = MakePatch();
  pml = MakePML({{"Attributes", {7}}});
  ExtrudeEntry(mesh, pml, {});
  CHECK(pml.extrusion.center ==
        std::array<std::array<double, 3>, 2>{{{0.0, 0.0, 0.0}, {0.0, 1.0, 0.5}}});
  mesh = MakePatch();
  pml = MakePML({{"Attributes", {7}}, {"Center", {0.0, 1.0, 0.5}}});
  CHECK_THROWS_WITH(ExtrudeEntry(mesh, pml, {}),
                    Matchers::ContainsSubstring("is not on its open side"));
}

TEST_CASE("Extruded PML depth", "[pml][geodata][Serial]")
{
  // Each form of the depth gives the depth of the nodes of the layer, from the position of
  // the surface on their direction.
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  const int order = GENERATE(1, 2);
  using Depth = config::PMLData::Extrusion::Depth;

  // The face z = 1 of a box, along its normal.
  auto smesh = std::make_unique<mfem::Mesh>(MakeBoxTop(type, 2, order));
  auto pml = MakePML({{"Attributes", {2}}});
  ExtrudeEntry(smesh, pml, {});
  REQUIRE(pml.extrusion.depth == Depth::NORMAL);
  CHECK(DepthError(*smesh, pml, [](const mfem::Vector &x) { return x(2) - 1.0; }) <
        1.0e-12);
  {
    // The uniformly refined layer keeps the depth, from 0 at the surface to the thickness
    // at the back.
    mfem::ParMesh pmesh(Mpi::World(), *smesh);
    pmesh.UniformRefinement();
    CHECK(DepthError(
              pmesh, pml, [](const mfem::Vector &x) { return x(2) - 1.0; }, false) <
          1.0e-12);
    double lo = mfem::infinity(), hi = -mfem::infinity();
    DepthError(
        pmesh, pml,
        [&](const mfem::Vector &x)
        {
          lo = std::min(lo, x(2) - 1.0);
          hi = std::max(hi, x(2) - 1.0);
          return x(2) - 1.0;
        },
        false);
    CHECK_THAT(lo, Matchers::WithinAbs(0.0, 1.0e-12));
    CHECK_THAT(hi, Matchers::WithinAbs(0.5, 1.0e-12));
  }

  // The same face of the box tilted by a rotation, along its tilted normal n, with the
  // depth n · x - 1, also translated by L (1, 0.7, 0.3) far from the origin, where n · x is
  // up to 10⁵ times the size of the mesh, with the depth n · x - 1 - L n · (1, 0.7, 0.3),
  // within the rounding of the coordinates.
  for (const double L : {0.0, 1.0e3, 1.0e5})
  {
    INFO("L = " << L);
    const double tol = 1.0e-12 + 1.0e-15 * L;
    smesh = std::make_unique<mfem::Mesh>(MakeBoxTop(type, 2, order));
    const double c = std::cos(0.4), s = std::sin(0.4);  // About the axis (1, 1, 1) / √3
    const double k = (1.0 - c) / 3.0, r = s / std::sqrt(3.0);
    const double R[3][3] = {
        {c + k, k - r, k + r}, {k + r, c + k, k - r}, {k - r, k + r, c + k}};
    const double shift[3] = {L, 0.7 * L, 0.3 * L};
    smesh->Transform(
        [&](const mfem::Vector &x, mfem::Vector &y)
        {
          y.SetSize(3);
          for (int i = 0; i < 3; i++)
          {
            y(i) = R[i][0] * x(0) + R[i][1] * x(1) + R[i][2] * x(2) + shift[i];
          }
        });
    const int ne = smesh->GetNE();
    pml = MakePML({{"Attributes", {2}}});
    ExtrudeEntry(smesh, pml, {});
    REQUIRE(pml.extrusion.depth == Depth::NORMAL);
    for (int i = 0; i < 3; i++)
    {
      CHECK_THAT(pml.extrusion.plane[i], Matchers::WithinAbs(R[i][2], tol));
    }
    const double h = 1.0 + R[0][2] * shift[0] + R[1][2] * shift[1] + R[2][2] * shift[2];
    CHECK(DepthError(*smesh, pml, [&](const mfem::Vector &x)
                     { return R[0][2] * x(0) + R[1][2] * x(1) + R[2][2] * x(2) - h; }) <
          tol);
    CHECK(MinExtrudedJacobian(*smesh, ne, order) > 0.0);
  }

  // A sphere, at the distance 1 from its centre.
  if (order == 2)
  {
    smesh = std::make_unique<mfem::Mesh>(MakeBall(type, 2, order));
    pml = MakePML();
    ExtrudeEntry(smesh, pml, {});
    REQUIRE(pml.extrusion.depth == Depth::DISTANCE);
    CHECK_THAT(pml.extrusion.distance, Matchers::WithinAbs(1.0, 1.0e-12));
    CHECK(DepthError(*smesh, pml, [](const mfem::Vector &x) { return x.Norml2() - 1.0; }) <
          1.0e-12);

    // The wall of a cylinder, at the distance 1 from its axis.
    smesh = std::make_unique<mfem::Mesh>(MakeCylinder(type, 2, order));
    pml = MakePML({{"Attributes", {1}}, {"Center", {{0.0, 0.0, -1.0}, {0.0, 0.0, 1.0}}}});
    ExtrudeEntry(smesh, pml, {});
    REQUIRE(pml.extrusion.depth == Depth::DISTANCE);
    CHECK(DepthError(*smesh, pml, [](const mfem::Vector &x)
                     { return std::hypot(x(0), x(1)) - 1.0; }) < 1.0e-12);
  }

  // A cube, extruded per component from its centre into a larger cube, and a layered box
  // from a segment through it, as Cartesian layers on all of their sides.
  smesh = std::make_unique<mfem::Mesh>(MakeBox(type, 2, order));
  pml = MakePML();
  ExtrudeEntry(smesh, pml, {});
  REQUIRE(pml.extrusion.depth == Depth::BOX);
  CHECK(std::ranges::all_of(pml.extrusion.box_sides, [](bool b) { return b; }));
  CHECK(BoxDepthError(*smesh, pml) < 1.0e-12);
  smesh = std::make_unique<mfem::Mesh>(MakeBox(type, 4, order, true));
  pml = MakePML({{"Center", {{0.0, 0.0, -0.5}, {0.0, 0.0, 0.5}}}});
  ExtrudeEntry(smesh, pml, {});
  REQUIRE(pml.extrusion.depth == Depth::BOX);
  CHECK(BoxDepthError(*smesh, pml) < 1.0e-12);
}

namespace
{

// Extrude the entries of config["Boundaries"]["PML"], replaced by the blocks of their
// parts.
std::vector<config::PMLData> ExtrudeEntries(std::unique_ptr<mfem::Mesh> &mesh,
                                            const std::vector<json> &entries)
{
  config::DomainData domains;
  for (std::size_t i = 0; i < entries.size(); i++)
  {
    domains.pml.emplace_back(entries[i], true).extrusion.entry = static_cast<int>(i);
  }
  mesh::ExtrudePML(mesh, domains, {});
  return domains.pml;
}

// Number of boundary elements with the given attribute, and whether every exterior face of
// the mesh has a boundary element.
int CountBdrElements(const mfem::Mesh &mesh, int attr)
{
  int count = 0;
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    count += (mesh.GetBdrAttribute(be) == attr);
  }
  return count;
}

bool ExteriorHasBdrElements(const mfem::Mesh &mesh)
{
  std::vector<bool> has_bdr(mesh.GetNumFaces(), false);
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    has_bdr[mesh.GetBdrElementFaceIndex(be)] = true;
  }
  for (int f = 0; f < mesh.GetNumFaces(); f++)
  {
    int e1, e2;
    mesh.GetFaceElements(f, &e1, &e2);
    if (e2 < 0 && !has_bdr[f])
    {
      return false;
    }
  }
  return true;
}

}  // namespace

TEST_CASE("ExtrudePML splits the surface into parts", "[pml][geodata][Serial]")
{
  using Depth = config::PMLData::Extrusion::Depth;
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  const int n = 2;

  SECTION("Pipe")
  {
    // The ends z = ±1 of the box (attributes 5 and 6) of one entry are two flat parts, each
    // along its own normal, or of two entries with their own thickness. The walls
    // (attribute 1) extend through both layers.
    for (const bool two : {false, true})
    {
      auto mesh = std::make_unique<mfem::Mesh>(MakeBox(type, n, 1));
      for (int be = 0; be < mesh->GetNBE(); be++)
      {
        mfem::Vector c(3);
        mfem::IsoparametricTransformation T;
        mesh->GetBdrElementTransformation(be, &T);
        T.Transform(mfem::Geometries.GetCenter(T.GetGeometryType()), c);
        mesh->SetBdrAttribute(be, (c(2) > 1.0 - 1.0e-12)    ? 6
                                  : (c(2) < -1.0 + 1.0e-12) ? 5
                                                            : 1);
      }
      mesh->SetAttributes();
      const int ne = mesh->GetNE(), n_walls = CountBdrElements(*mesh, 1);
      const json end = {{"Thickness", 0.1}, {"Layers", 2}};
      auto blocks =
          two ? ExtrudeEntries(mesh,
                               {{{"Attributes", {5}}, {"Thickness", 0.1}, {"Layers", 2}},
                                {{"Attributes", {6}}, {"Thickness", 0.2}, {"Layers", 3}}})
              : ExtrudeEntries(
                    mesh, {{{"Attributes", {5, 6}}, {"Thickness", 0.1}, {"Layers", 2}}});
      REQUIRE(blocks.size() == 2);
      for (const auto &block : blocks)
      {
        REQUIRE(block.extrusion.depth == Depth::NORMAL);
        const double nz = block.extrusion.plane[2];
        CHECK_THAT(std::abs(nz), Matchers::WithinAbs(1.0, 1.0e-12));
        CHECK(DepthError(*mesh, block,
                         [&](const mfem::Vector &x) { return nz * x(2) - 1.0; }) < 1.0e-12);
      }
      CHECK(blocks[0].attributes != blocks[1].attributes);
      CHECK(blocks[0].extrusion.back_attribute == blocks[1].extrusion.back_attribute);
      CHECK(MinExtrudedJacobian(*mesh, ne, 1) > 0.0);
      // 4 n edges at each end, each extended into one quadrilateral per layer.
      CHECK(CountBdrElements(*mesh, 1) == n_walls + 4 * n * (two ? 2 + 3 : 2 * 2));
      CHECK(ExteriorHasBdrElements(*mesh));
    }
  }

  SECTION("Cube")
  {
    // By default, or with the radial extrusion, the cube is one part extruded per
    // component from its centre, or from a box just inside of it (straight over the box,
    // with fans at the edges and corners), into a larger cube: the layer has the thickness
    // normal to every face. The normal extrusion, for separate flat faces, would make each
    // face a part, which shares its edges with four others.
    for (const std::string extrusion : {"", "Radial", "Box"})
    {
      auto mesh = std::make_unique<mfem::Mesh>(MakeBox(type, 3, 1));
      const int ne = mesh->GetNE(), nbe = mesh->GetNBE();
      json entry = {{"Attributes", {1}}, {"Thickness", 0.5}, {"Layers", 2}};
      if (extrusion == "Box")
      {
        entry["Center"] = {{-0.9, -0.9, -0.9}, {0.9, 0.9, 0.9}};
      }
      else if (!extrusion.empty())
      {
        entry["Extrusion"] = extrusion;
      }
      auto blocks = ExtrudeEntries(mesh, {entry});
      CHECK(MinExtrudedJacobian(*mesh, ne, 1) > 0.0);
      CHECK(ExteriorHasBdrElements(*mesh));
      const int back = blocks[0].extrusion.back_attribute;
      REQUIRE(blocks.size() == 1);
      CHECK(blocks[0].extrusion.depth == Depth::BOX);
      CHECK(CountBdrElements(*mesh, back) == nbe);
      CHECK(BoxDepthError(*mesh, blocks[0]) < 1.0e-12);
      double err = 0.0;
      for (int be = nbe; be < mesh->GetNBE(); be++)
      {
        for (int j = 0; mesh->GetBdrAttribute(be) == back &&
                        j < mesh->GetNodes()->FESpace()->GetBE(be)->GetDof();
             j++)
        {
          err = std::max(err, std::abs(NodePosition(*mesh, be, j, true).Normlinf() - 1.5));
        }
      }
      CHECK(err < 1.0e-12);
    }
    auto mesh = std::make_unique<mfem::Mesh>(MakeBox(type, 3, 1));
    CHECK_THROWS_WITH(
        ExtrudeEntries(
            mesh, {{{"Attributes", {1}}, {"Thickness", 0.5}, {"Extrusion", "Normal"}}}),
        Matchers::ContainsSubstring("which share the edge at") &&
            Matchers::ContainsSubstring("is for separate flat faces"));
  }

  SECTION("Cylinder")
  {
    // The wall (attribute 1) is at the distance 1 from its axis, and its layer ends flush
    // with the caps (attribute 2), which continue through it.
    auto mesh = std::make_unique<mfem::Mesh>(MakeCylinder(type, n, 2));
    const int ne = mesh->GetNE(), n_caps = CountBdrElements(*mesh, 2);
    auto blocks = ExtrudeEntries(mesh, {{{"Attributes", {1}},
                                         {"Thickness", 0.1},
                                         {"Layers", 2},
                                         {"Center", {{0.0, 0.0, -1.0}, {0.0, 0.0, 1.0}}}}});
    REQUIRE(blocks.size() == 1);
    CHECK(MinExtrudedJacobian(*mesh, ne, 2) > 0.0);
    CHECK(ExteriorHasBdrElements(*mesh));
    REQUIRE(blocks[0].extrusion.depth == Depth::DISTANCE);
    CHECK(DepthError(*mesh, blocks[0], [](const mfem::Vector &x)
                     { return std::hypot(x(0), x(1)) - 1.0; }) < 1.0e-12);
    // The rim of each cap has 4 n edges, with one quadrilateral per layer.
    CHECK(CountBdrElements(*mesh, 2) == n_caps + 2 * (4 * n * 2));

    // The closed cylinder, with its caps as parts of their own, whose layers would meet the
    // wall's along the rims.
    mesh = std::make_unique<mfem::Mesh>(MakeCylinder(type, n, 2));
    CHECK_THROWS_WITH(
        ExtrudeEntries(mesh, {{{"Attributes", {2}}, {"Thickness", 0.1}},
                              {{"Attributes", {1}},
                               {"Thickness", 0.1},
                               {"Center", {{0.0, 0.0, -1.0}, {0.0, 0.0, 1.0}}}}}),
        Matchers::ContainsSubstring("which share the edge at") &&
            Matchers::ContainsSubstring("list their attributes in one entry, so that"));

    // In one entry, the closed cylinder, of flat and curved faces, fits no shape with any
    // core: the error says so, and what to extrude instead.
    for (const auto &center : {json(), json({{0.0, 0.0, -1.0}, {0.0, 0.0, 1.0}})})
    {
      mesh = std::make_unique<mfem::Mesh>(MakeCylinder(type, n, 2));
      json entry = {{"Attributes", {1, 2}}, {"Thickness", 0.1}};
      if (!center.is_null())
      {
        entry["Center"] = center;
      }
      CHECK_THROWS_WITH(
          ExtrudeEntries(mesh, {entry}),
          Matchers::ContainsSubstring("which has flat and curved faces") &&
              Matchers::ContainsSubstring("for example a closed cylinder, cannot be") &&
              Matchers::ContainsSubstring("from a box or a sphere around the object"));
    }

    // The wall is curved, and cannot be extruded along its normal.
    mesh = std::make_unique<mfem::Mesh>(MakeCylinder(type, n, 2));
    CHECK_THROWS_WITH(
        ExtrudeEntries(
            mesh, {{{"Attributes", {1}}, {"Thickness", 0.1}, {"Extrusion", "Normal"}}}),
        Matchers::ContainsSubstring("is curved"));
  }

  SECTION("Two entries on one flat face")
  {
    // The face z = 1 of the box split at x = 0 into the attributes 3 and 5 is one flat part
    // of one entry, and two parts which share the edges x = 0 of two entries, whose layers
    // would be split by their sides.
    auto MakeSplitTop = [&]()
    {
      auto mesh = std::make_unique<mfem::Mesh>(MakeBoxTop(type, n, 1));
      for (int be = 0; be < mesh->GetNBE(); be++)
      {
        mfem::Vector c(3);
        mfem::IsoparametricTransformation T;
        mesh->GetBdrElementTransformation(be, &T);
        T.Transform(mfem::Geometries.GetCenter(T.GetGeometryType()), c);
        if (mesh->GetBdrAttribute(be) == 2)
        {
          mesh->SetBdrAttribute(be, (c(0) < 0.0) ? 3 : 5);
        }
      }
      mesh->SetAttributes();
      return mesh;
    };
    auto mesh = MakeSplitTop();
    const auto blocks =
        ExtrudeEntries(mesh, {{{"Attributes", {3, 5}}, {"Thickness", 0.1}, {"Layers", 2}}});
    REQUIRE(blocks.size() == 1);
    CHECK(blocks[0].extrusion.depth == Depth::NORMAL);
    CHECK(ExteriorHasBdrElements(*mesh));
    mesh = MakeSplitTop();
    CHECK_THROWS_WITH(
        ExtrudeEntries(mesh, {{{"Attributes", {3}}, {"Thickness", 0.1}},
                              {{"Attributes", {5}}, {"Thickness", 0.1}}}),
        Matchers::ContainsSubstring("which share the edge at (0, ") &&
            Matchers::ContainsSubstring("list their attributes in one entry"));
  }

  SECTION("Two parts which meet at a vertex")
  {
    // The quadrilaterals x, y > 0 of the face z = 1 and y < 0, z > 0 of the face x = 1
    // (attributes 3 and 4) meet at the vertex (1, 0, 1) only.
    auto mesh = std::make_unique<mfem::Mesh>(MakeBox(mfem::Element::HEXAHEDRON, n, 1));
    for (int be = 0; be < mesh->GetNBE(); be++)
    {
      mfem::Vector c(3);
      mfem::IsoparametricTransformation T;
      mesh->GetBdrElementTransformation(be, &T);
      T.Transform(mfem::Geometries.GetCenter(T.GetGeometryType()), c);
      mesh->SetBdrAttribute(be, (c(2) > 1.0 - 1.0e-12 && c(0) > 0.0 && c(1) > 0.0)   ? 3
                                : (c(0) > 1.0 - 1.0e-12 && c(1) < 0.0 && c(2) > 0.0) ? 4
                                                                                     : 1);
    }
    mesh->SetAttributes();
    CHECK_THROWS_WITH(
        ExtrudeEntries(mesh, {{{"Attributes", {3, 4}}, {"Thickness", 0.5}, {"Layers", 2}}}),
        Matchers::ContainsSubstring(
            "which meet at the vertex (1, 0, 1)! Two extruded parts "
            "must not touch"));
  }
}

TEST_CASE("ExtrudePML along the normal of a periodic boundary", "[pml][geodata][Serial]")
{
  // The top z = 1 of the cell [0, 1]³, periodic in x and y: the layer along its normal is
  // periodic too, with each periodic boundary extended through it on the faces between its
  // columns across the seams, as two boundary elements of the donor and receiver attributes
  // on each face.
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  const int order = GENERATE(1, 2);
  auto cell = mfem::Mesh::MakeCartesian3D(3, 3, 2, type, 1.0, 1.0, 1.0);
  cell.SetCurvature(order);
  std::vector<mfem::Vector> translations(2, mfem::Vector(3));
  for (int i = 0; i < 2; i++)
  {
    translations[i] = 0.0;
    translations[i](i) = 1.0;
  }
  const auto periodic_mesh =
      mfem::Mesh::MakePeriodic(cell, cell.CreatePeriodicVertexMapping(translations));
  auto Extrude = [&](json extra)
  {
    // Boundary attributes of mfem::Mesh::MakeCartesian3D: 1 to 6 for z = 0, y = 0, x = 1,
    // y = 1, x = 0, z = 1.
    config::BoundaryData boundaries;
    boundaries.periodic.boundary_pairs.resize(2);
    boundaries.periodic.boundary_pairs[0].donor_attributes = {5};
    boundaries.periodic.boundary_pairs[0].receiver_attributes = {3};
    boundaries.periodic.boundary_pairs[1].donor_attributes = {2};
    boundaries.periodic.boundary_pairs[1].receiver_attributes = {4};
    auto mesh = std::make_unique<mfem::Mesh>(periodic_mesh);
    auto pml = MakePML(extra);
    ExtrudeEntry(mesh, pml, boundaries);
    return std::make_pair(std::move(mesh), pml);
  };
  const int ne = periodic_mesh.GetNE(), nbe = periodic_mesh.GetNBE();
  auto [mesh, pml] = Extrude({{"Attributes", {6}}});
  REQUIRE(pml.extrusion.depth == config::PMLData::Extrusion::Depth::NORMAL);
  CHECK(DepthError(*mesh, pml, [](const mfem::Vector &x) { return x(2) - 1.0; }) < 1.0e-12);
  CHECK(MinExtrudedJacobian(*mesh, ne, order) > 0.0);
  std::map<int, int> count;
  for (int be = nbe; be < mesh->GetNBE(); be++)
  {
    const int attr = mesh->GetBdrAttribute(be);
    if (attr == pml.extrusion.back_attribute)
    {
      continue;
    }
    count[attr]++;
    int e1, e2;
    mesh->GetFaceElements(mesh->GetBdrElementFaceIndex(be), &e1, &e2);
    CHECK(e1 >= ne);
    CHECK(e2 >= ne);
  }
  CHECK(count == std::map<int, int>{{2, 6}, {3, 6}, {4, 6}, {5, 6}});

  // Along rays, the two sides of the periodic boundaries would flare differently, as for a
  // top curved into a cylinder of radius 2 around an axis along y.
  {
    auto curved = mfem::Mesh::MakeCartesian3D(3, 3, 2, type, 1.0, 1.0, 1.0);
    curved.SetCurvature(2);
    curved.Transform(
        [](const mfem::Vector &x, mfem::Vector &y)
        {
          y = x;
          y(2) = x(2) * (std::sqrt(4.0 - (x(0) - 0.5) * (x(0) - 0.5)) - 1.0);
        });
    auto mesh = std::make_unique<mfem::Mesh>(
        mfem::Mesh::MakePeriodic(curved, curved.CreatePeriodicVertexMapping(translations)));
    config::BoundaryData boundaries;
    boundaries.periodic.boundary_pairs.resize(1);
    boundaries.periodic.boundary_pairs[0].donor_attributes = {5};
    boundaries.periodic.boundary_pairs[0].receiver_attributes = {3};
    auto pml =
        MakePML({{"Attributes", {6}}, {"Center", {{0.5, 0.0, -1.0}, {0.5, 1.0, -1.0}}}});
    config::DomainData domains;
    domains.pml = {pml};
    CHECK_THROWS_WITH(mesh::ExtrudePML(mesh, domains, boundaries),
                      Matchers::ContainsSubstring("would flare the two sides"));
  }
}

namespace
{

// Extruded stretch at the distance 1 from a point core at the origin, with the profile of a
// layer of thickness 0.5.
pml::Stretch MakeExtrudedStretch(bool frequency_dependent)
{
  pml::Stretch p;
  p.extruded = true;
  p.depth = config::PMLData::Extrusion::Depth::DISTANCE;
  p.distance = 1.0;
  p.geometry.thickness[0] = 0.5;
  p.sigma_max[0] = 3.0;
  p.kappa_max = {2.0, 2.0, 2.0};
  p.alpha_max.fill(frequency_dependent ? 0.4 : 0.0);
  p.order = 3;
  p.frequency_dependent = frequency_dependent;
  p.reference_frequency = frequency_dependent ? 0.0 : 1.5;
  return p;
}

// Independent evaluation of the stretch s(ξ) and of ξ̃(ξ) = ∫₀^ξ s (by Simpson's rule).
std::pair<std::complex<double>, std::complex<double>>
ReferenceProfile(const pml::Stretch &p, double xi, std::complex<double> omega)
{
  if (!p.frequency_dependent)
  {
    omega = p.reference_frequency;
  }
  // Beyond the layer, the stretch of its outer edge.
  const double t = p.geometry.thickness[0];
  auto s = [&](double x)
  {
    const double shape = std::pow(std::min(x, t) / t, p.order);
    return 1.0 + (p.kappa_max[0] - 1.0) * shape +
           p.sigma_max[0] * shape / (p.alpha_max[0] * shape + 1i * omega);
  };
  constexpr int n = 2000;
  const double xi_t = std::min(xi, t);
  std::complex<double> integral = (xi - xi_t) * s(t);
  for (int i = 0; i < n; i++)
  {
    const double a = xi_t * i / n, b = xi_t * (i + 1) / n;
    integral += (b - a) / 6.0 * (s(a) + 4.0 * s(0.5 * (a + b)) + s(b));
  }
  return {s(xi), integral};
}

// Requested part of c μ̃⁻¹ or c ε̃ from the device helpers, for an isotropic or anisotropic
// background.
std::array<double, 9> EvalCoeffAt(const pml::Stretch &p, const pml::Background &b,
                                  pml::TensorPart part, std::complex<double> c,
                                  std::complex<double> omega, bool muinv,
                                  const std::array<double, 3> &x)
{
  pml::ContextHeader header;
  header.part = part;
  (muinv ? header.c_muinv : header.c_eps) = c;
  header.omega = omega;
  const auto ctx = pml::PackContext(header, p, {0}, {b});
  std::array<double, 9> coeff;
  REQUIRE(PMLCoeffAt(ctx.data(), 1, x.data(), muinv ? coeff.data() : nullptr,
                     muinv ? nullptr : coeff.data()));
  return coeff;
}

double ExpectedPart(pml::TensorPart part, std::complex<double> c, std::complex<double> t,
                    std::complex<double> t_mag)
{
  switch (part)
  {
    case pml::TensorPart::REAL:
      return (c * t).real();
    case pml::TensorPart::IMAG:
      return (c * t).imag();
    case pml::TensorPart::STRETCH_MAGNITUDE_REAL:
      return (c * t_mag).real();
    case pml::TensorPart::STRETCH_MAGNITUDE_IMAG:
      return (c * t_mag).imag();
  }
  return 0.0;
}

constexpr std::array<pml::TensorPart, 4> all_parts = {
    pml::TensorPart::REAL, pml::TensorPart::IMAG, pml::TensorPart::STRETCH_MAGNITUDE_REAL,
    pml::TensorPart::STRETCH_MAGNITUDE_IMAG};

}  // namespace

TEST_CASE("Extruded PML stretch", "[pml][Serial]")
{
  const bool frequency_dependent = GENERATE(false, true);
  const auto p = MakeExtrudedStretch(frequency_dependent);
  const std::complex<double> omega = 1.3 + 0.05i, c = 0.7 - 0.2i;

  SECTION("Identity at the interface")
  {
    // On the surface, and inside of it (a negative depth, of a point of a sampled curved
    // surface, is zero).
    for (double r : {1.0, 0.98})
    {
      const std::array<double, 3> x = {0.48 * r, -0.6 * r, 0.64 * r};
      for (bool magnitude : {false, true})
      {
        const auto J = pml::EvaluateJacobian(p, x, omega, magnitude);
        for (int k = 0; k < 9; k++)
        {
          CHECK_THAT(J[k].real(), Matchers::WithinAbs((k % 4 == 0) ? 1.0 : 0.0, 1.0e-15));
          CHECK(J[k].imag() == 0.0);
        }
      }
    }
  }

  SECTION("Spherical PML")
  {
    // On the sphere of radius 1 around the core, the depth is ξ = r - 1 with ∇ξ = d, and
    // the tensors are those of the spherical PML for the radial stretch s and the
    // tangential stretch r̃ / r, r̃ = 1 + ξ̃: μ̃⁻¹ = μ⁻¹ (s / (r̃ / r)² d dᵀ + (1 / s) (I - d
    // dᵀ)) and ε̃ = ε ((r̃ / r)² / s d dᵀ + s (I - d dᵀ)).
    pml::Background b;
    const double mu_inv = 1.0 / 1.5;
    const std::complex<double> eps = 2.5 - 0.01i;
    for (int i = 0; i < 3; i++)
    {
      b.mu_inv[4 * i] = mu_inv;
      b.epsilon_real[4 * i] = eps.real();
      b.epsilon_imag[4 * i] = eps.imag();
    }
    const std::array<double, 3> dir = {0.48, -0.6, 0.64};
    for (double xi : {0.1, 0.25, 0.4, 0.5})
    {
      const double r = 1.0 + xi;
      const std::array<double, 3> x = {r * dir[0], r * dir[1], r * dir[2]};
      const auto [s, xi_t] = ReferenceProfile(p, xi, omega);
      const auto q = (1.0 + xi_t) / r;
      const double s_mag = std::abs(s), q_mag = std::abs(q);
      const auto s_x = pml::EvaluateStretch(p, x, omega);
      CHECK(std::abs(s_x[0] - s) < 1.0e-12);
      CHECK(std::abs(s_x[1] - q) < 1.0e-9);
      CHECK(std::abs(s_x[2] - q) < 1.0e-9);
      for (bool muinv : {true, false})
      {
        for (auto part : all_parts)
        {
          const auto coeff = EvalCoeffAt(p, b, part, c, omega, muinv, x);
          for (int j = 0; j < 3; j++)
          {
            for (int i = 0; i < 3; i++)
            {
              const double dd = dir[i] * dir[j], tt = (i == j) - dd;
              const auto t = muinv ? mu_inv * (s / (q * q) * dd + tt / s)
                                   : eps * (q * q / s * dd + s * tt);
              const auto t_mag = muinv
                                     ? mu_inv * (s_mag / (q_mag * q_mag) * dd + tt / s_mag)
                                     : eps * (q_mag * q_mag / s_mag * dd + s_mag * tt);
              INFO("xi = " << xi << ", muinv = " << muinv
                           << ", part = " << static_cast<int>(part));
              CHECK_THAT(coeff[i + 3 * j],
                         Matchers::WithinAbs(ExpectedPart(part, c, t, t_mag), 1.0e-9));
            }
          }
        }
      }
    }
  }

  SECTION("Tangential stretch with a large frequency shift")
  {
    // J̃ = s d dᵀ + (r̃ / r) (I - d dᵀ) at the outer edge of the layer for α / ω = 1, 3 and
    // 10 (σ / ω = 20, κ = 1, real ω), where ξ̃ by 8-point Gauss-Legendre quadrature was off
    // by 8e-7, 5e-5 and 3e-4 in r̃ / r.
    auto p_a = MakeExtrudedStretch(true);
    p_a.sigma_max[0] = 20.0;
    p_a.kappa_max.fill(1.0);
    const std::array<double, 3> dir = {0.48, -0.6, 0.64};
    const double xi = p_a.geometry.thickness[0], r = 1.0 + xi;
    const std::array<double, 3> x = {r * dir[0], r * dir[1], r * dir[2]};
    for (double alpha : {1.0, 3.0, 10.0})
    {
      p_a.alpha_max.fill(alpha);
      const auto [s, xi_t] = ReferenceProfile(p_a, xi, 1.0);
      const auto q = (1.0 + xi_t) / r;
      const auto J = pml::EvaluateJacobian(p_a, x, 1.0, false);
      for (int j = 0; j < 3; j++)
      {
        for (int i = 0; i < 3; i++)
        {
          const double dd = dir[i] * dir[j];
          INFO("alpha = " << alpha << ", i = " << i << ", j = " << j);
          CHECK(std::abs(J[i + 3 * j] - (s * dd + q * ((i == j) - dd))) < 1.0e-12);
        }
      }
    }
  }

  SECTION("Cartesian PML on a face of a box")
  {
    // For normal extrusion of a face (d ∇ξᵀ = n nᵀ, ∇d = 0), J̃ = I + (s - 1) n nᵀ is the
    // diagonal stretch S of the Cartesian PML of the face, and the tensors are the same for
    // an anisotropic, lossy background.
    auto p_n = p;
    p_n.depth = config::PMLData::Extrusion::Depth::NORMAL;
    p_n.plane = {0.0, 0.0, 1.0, 1.0};
    pml::Stretch cartesian;
    cartesian.geometry.inner = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
    cartesian.geometry.thickness = {0.0, 0.0, 0.0, 0.0, 0.0, p.geometry.thickness[0]};
    cartesian.sigma_max[5] = p.sigma_max[0];
    cartesian.kappa_max = p.kappa_max;
    cartesian.alpha_max = p.alpha_max;
    cartesian.order = p.order;
    cartesian.frequency_dependent = p.frequency_dependent;
    cartesian.reference_frequency = p.reference_frequency;
    pml::Background b;
    b.mu_inv = {0.8, 0.1, 0.0, 0.1, 0.9, 0.05, 0.0, 0.05, 1.1};
    b.epsilon_real = {4.0, 0.3, 0.2, 0.3, 5.0, 0.1, 0.2, 0.1, 6.0};
    b.epsilon_imag = {-4e-3, 0.0, 0.0, 0.0, -5e-3, 0.0, 0.0, 0.0, -6e-3};
    for (double xi : {0.1, 0.3, 0.5})
    {
      const std::array<double, 3> x = {0.2, -0.4, 1.0 + xi};
      const auto s = pml::EvaluateStretch(cartesian, x, omega)[2];
      const auto s_n = pml::EvaluateStretch(p_n, x, omega);
      CHECK(std::abs(s_n[0] - s) < 1.0e-12);
      CHECK(s_n[1] == 1.0);
      CHECK(s_n[2] == 1.0);
      for (bool muinv : {true, false})
      {
        for (auto part : all_parts)
        {
          pml::ContextHeader header;
          header.part = part;
          (muinv ? header.c_muinv : header.c_eps) = c;
          header.omega = omega;
          const auto ctx_s = pml::PackContext(header, cartesian, {0}, {b});
          const auto ctx_n = pml::PackContext(header, p_n, {0}, {b});
          std::array<double, 9> coeff_s, coeff_j;
          PMLCoeffAt(ctx_s.data(), 1, x.data(), muinv ? coeff_s.data() : nullptr,
                     muinv ? nullptr : coeff_s.data());
          PMLCoeffAt(ctx_n.data(), 1, x.data(), muinv ? coeff_j.data() : nullptr,
                     muinv ? nullptr : coeff_j.data());
          REQUIRE(std::abs(s - 1.0) > 0.01);
          for (int k = 0; k < 9; k++)
          {
            CHECK_THAT(coeff_j[k], Matchers::WithinAbs(coeff_s[k], 1.0e-12));
          }
        }
      }
    }
  }
}

TEST_CASE("Extruded PML depth gradients", "[pml][Serial]")
{
  // ∇ξ, ∇d, and J̃ of each form of the depth against central differences, at points of the
  // layer away from the kinks of the depth (the ends of the core), for the static and the
  // frequency-dependent stretch (with α > 0, whose ξ̃ is by quadrature). The points deeper
  // in the layer (ξ = 0.35) are where ξ̃ - ξ, and its α term, weigh in J̃.
  using Depth = config::PMLData::Extrusion::Depth;
  struct Case
  {
    Depth depth;
    std::array<std::array<double, 3>, 2> center;
    std::array<double, 3> x;
  };
  // Plain variables, since Clang 18 cannot capture a structured binding with OpenMP.
  const auto test = GENERATE_COPY(
      Case{Depth::NORMAL, {}, {0.3, -0.2, 1.2}},
      Case{Depth::DISTANCE, {{{0.1, 0.0, 0.0}, {0.1, 0.0, 0.0}}}, {0.5, -0.9, 0.4}},
      Case{Depth::DISTANCE, {{{0.0, 0.0, -1.0}, {0.0, 0.0, 1.0}}}, {0.6, -0.9, 0.4}},
      Case{Depth::DISTANCE, {{{-0.5, -0.2, 0.0}, {0.5, 0.2, 0.0}}}, {0.9, -1.0, 0.6}},
      Case{Depth::DISTANCE, {{{0.1, 0.0, 0.0}, {0.1, 0.0, 0.0}}}, {0.748, -0.81, 0.864}},
      Case{Depth::DISTANCE, {{{0.0, 0.0, -1.0}, {0.0, 0.0, 1.0}}}, {0.749, -1.123, 0.4}},
      Case{
          Depth::DISTANCE, {{{-0.5, -0.2, 0.0}, {0.5, 0.2, 0.0}}}, {1.001, -1.203, 0.752}});
  const auto &x = test.x;
  const bool frequency_dependent = GENERATE(false, true);
  const std::complex<double> omega = 1.3 + 0.05i;
  auto p = MakeExtrudedStretch(frequency_dependent);
  p.depth = test.depth;
  p.center = test.center;
  p.distance = 1.0;
  p.plane = {0.0, 0.0, 1.0, 1.0};
  auto X = [&](int c, double h)
  {
    auto y = x;
    y[c] += h;
    return y;
  };
  // x̃ = x + (ξ̃(ξ) - ξ) d.
  auto Stretched = [&](const std::array<double, 3> &y)
  {
    const auto dy = pml::EvaluateDepth(p, y);
    const auto xi_t = ReferenceProfile(p, dy.xi, omega).second;
    std::array<std::complex<double>, 3> y_t;
    for (int c = 0; c < 3; c++)
    {
      y_t[c] = y[c] + (xi_t - dy.xi) * dy.d[c];
    }
    return y_t;
  };
  const auto d0 = pml::EvaluateDepth(p, x);
  REQUIRE(d0.xi > 0.05);
  REQUIRE(d0.xi < 0.45);
  const auto J = pml::EvaluateJacobian(p, x, omega);
  constexpr double h = 1.0e-5;
  for (int j = 0; j < 3; j++)
  {
    const auto dp = pml::EvaluateDepth(p, X(j, h)), dm = pml::EvaluateDepth(p, X(j, -h));
    CHECK_THAT(d0.grad_xi[j], Matchers::WithinAbs((dp.xi - dm.xi) / (2.0 * h), 1.0e-8));
    const auto yp = Stretched(X(j, h)), ym = Stretched(X(j, -h));
    for (int i = 0; i < 3; i++)
    {
      CHECK_THAT(d0.grad_d[i + 3 * j],
                 Matchers::WithinAbs((dp.d[i] - dm.d[i]) / (2.0 * h), 1.0e-8));
      const auto J_fd = (yp[i] - ym[i]) / (2.0 * h);
      CHECK_THAT(J[i + 3 * j].real(), Matchers::WithinAbs(J_fd.real(), 1.0e-6));
      CHECK_THAT(J[i + 3 * j].imag(), Matchers::WithinAbs(J_fd.imag(), 1.0e-6));
    }
  }
}

namespace
{

// The tensors of the extruded stretch at a quadrature point, evaluated independently of the
// device helpers, for a layer at the distance r₀ from a point core at the origin (the depth
// ξ = r - r₀, ∇ξ = x / r) and an isotropic background, with the frequency-dependent stretch
// evaluated at omega (ignored for a static stretch).
class ReferenceExtrudedTensor : public mfem::MatrixCoefficient
{
private:
  const pml::Stretch &p;
  double mu_inv;
  std::complex<double> eps, omega;
  bool muinv, imag;

public:
  ReferenceExtrudedTensor(const pml::Stretch &p, double mu_inv, std::complex<double> eps,
                          std::complex<double> omega, bool muinv, bool imag)
    : mfem::MatrixCoefficient(3), p(p), mu_inv(mu_inv), eps(eps), omega(omega),
      muinv(muinv), imag(imag)
  {
  }

  void Eval(mfem::DenseMatrix &K, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    T.SetIntPoint(&ip);
    mfem::Vector x(3), grad_xi(x);
    T.Transform(ip, x);
    const double r = x.Norml2(), xi = std::max(r - p.distance, 0.0);
    grad_xi = x;
    grad_xi /= r;
    const auto [s, xi_t] = ReferenceProfile(p, xi, omega);
    std::array<std::complex<double>, 9> J, adj;
    for (int j = 0; j < 3; j++)
    {
      for (int i = 0; i < 3; i++)
      {
        const double d_i = x(i) / r, grad_d = ((i == j) - x(i) * x(j) / (r * r)) / r;
        J[i + 3 * j] = double(i == j) + (s - 1.0) * d_i * grad_xi(j) + (xi_t - xi) * grad_d;
      }
    }
    for (int j = 0; j < 3; j++)
    {
      for (int i = 0; i < 3; i++)
      {
        const int i1 = (i + 1) % 3, i2 = (i + 2) % 3, j1 = (j + 1) % 3, j2 = (j + 2) % 3;
        adj[i + 3 * j] =
            J[j1 + 3 * i1] * J[j2 + 3 * i2] - J[j1 + 3 * i2] * J[j2 + 3 * i1];  // cofactor
      }
    }
    const auto det = J[0] * adj[0] + J[3] * adj[1] + J[6] * adj[2];
    for (int j = 0; j < 3; j++)
    {
      for (int i = 0; i < 3; i++)
      {
        std::complex<double> t = 0.0;
        for (int k = 0; k < 3; k++)
        {
          t += muinv ? J[k + 3 * i] * mu_inv * J[k + 3 * j]
                     : adj[i + 3 * k] * eps * adj[j + 3 * k];
        }
        t /= det;
        K(i, j) = imag ? t.imag() : t.real();
      }
    }
  }
};

}  // namespace

TEST_CASE("SpaceOperator extruded PML matches an MFEM coefficient reference",
          "[spaceoperator][pml][Serial][Parallel]")
{
  // Ball of hexahedra with an extruded spherical layer, also uniformly refined: the
  // libCEED operators of the PML terms match MFEM integrators with the tensors of the
  // spherical layer, for the static stretch and the frequency-dependent one with α > 0 at a
  // complex frequency.
  MPI_Comm comm = Mpi::World();
  const int refine = GENERATE(0, 1);
  const bool frequency_dependent = GENERATE(false, true);
  const std::complex<double> omega = 1.3 + 0.05i;
  auto serial_mesh =
      std::make_unique<mfem::Mesh>(MakeBall(mfem::Element::HEXAHEDRON, 2, 2));
  config::DomainData domains;
  domains.pml = {MakePML()};
  domains.pml[0].reference_frequency = 1.5;
  domains.pml[0].kappa_max.fill(1.5);
  domains.pml[0].frequency_dependent = frequency_dependent;
  domains.pml[0].alpha_max.fill(frequency_dependent ? 0.3 : 0.0);
  mesh::ExtrudePML(serial_mesh, domains, {});
  std::vector<std::unique_ptr<Mesh>> mesh;
  auto pmesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  for (int l = 0; l < refine; l++)
  {
    pmesh->UniformRefinement();
  }
  mesh.push_back(std::make_unique<Mesh>(std::move(pmesh)));

  const int pa_order_threshold = BilinearForm::pa_order_threshold;
  const int p_trial = fem::DefaultIntegrationOrder::p_trial;
  config::SolverData solver;
  solver.order = 2;
  solver.linear.mg_max_levels = 1;
  fem::DefaultIntegrationOrder::p_trial = solver.order;
  config::MaterialData material;
  material.attributes = {1, 2};
  material.mu_r.s.fill(1.5);
  material.epsilon_r.s.fill(2.5);
  domains.attributes = {1, 2};
  domains.materials = {material};
  config::BoundaryData boundaries;
  Units units(1.0, 1.0);
  SpaceOperator space_op(solver, domains, boundaries, ProblemType::DRIVEN, units, mesh);
  REQUIRE(space_op.GetMaterialOp().HasPML());
  const auto &p = space_op.GetMaterialOp().GetPMLLayers().at(0).GetStretch();
  REQUIRE(p.extruded);
  REQUIRE(p.frequency_dependent == frequency_dependent);

  const auto &pfes = space_op.GetNDSpace().Get();
  const auto &ir = mfem::IntRules.Get(
      mfem::Geometry::CUBE,
      fem::DefaultIntegrationOrder::Get(mesh[0]->Get(), mfem::Geometry::CUBE));
  const std::complex<double> eps = 2.5;
  auto AssembleReference = [&](bool muinv, bool imag)
  {
    mfem::ParBilinearForm a(const_cast<mfem::ParFiniteElementSpace *>(&pfes));
    mfem::Array<int> pml_marker({0, 1}), bulk_marker({1, 0});
    ReferenceExtrudedTensor pml_coef(p, 1.0 / 1.5, eps, omega, muinv, imag);
    mfem::ConstantCoefficient bulk_coef(muinv ? 1.0 / 1.5 : 2.5);
    auto *integ_pml = muinv ? static_cast<mfem::BilinearFormIntegrator *>(
                                  new mfem::CurlCurlIntegrator(pml_coef))
                            : new mfem::VectorFEMassIntegrator(pml_coef);
    integ_pml->SetIntRule(&ir);
    a.AddDomainIntegrator(integ_pml, pml_marker);
    if (!imag)
    {
      auto *integ_bulk = muinv ? static_cast<mfem::BilinearFormIntegrator *>(
                                     new mfem::CurlCurlIntegrator(bulk_coef))
                               : new mfem::VectorFEMassIntegrator(bulk_coef);
      integ_bulk->SetIntRule(&ir);
      a.AddDomainIntegrator(integ_bulk, bulk_marker);
    }
    a.Assemble();
    a.Finalize();
    return std::unique_ptr<mfem::HypreParMatrix>(a.ParallelAssemble());
  };

  // (Ar + i Ai) x for the reference of the μ̃⁻¹ or ε̃ terms.
  auto ApplyReference = [&](bool muinv, const ComplexVector &x, ComplexVector &y)
  {
    auto Ar = AssembleReference(muinv, false), Ai = AssembleReference(muinv, true);
    Vector t(y.Size());
    Ar->Mult(x.Real(), y.Real());
    Ai->Mult(x.Imag(), t);
    y.Real().Add(-1.0, t);
    Ar->Mult(x.Imag(), y.Imag());
    Ai->Mult(x.Real(), t);
    y.Imag().Add(1.0, t);
  };
  ComplexVector x(space_op.GetNDSpace().GetTrueVSize()), y(x.Size()), y_ref(x.Size()),
      t(x.Size());
  x.Real().Randomize(42 + Mpi::Rank(comm));
  x.Imag().Randomize(4242 + Mpi::Rank(comm));
  auto RelErr = [&]()
  {
    const double norm = linalg::Norml2(comm, y_ref);
    y.Add(-1.0, y_ref);
    return linalg::Norml2(comm, y) / norm;
  };
  if (!frequency_dependent)
  {
    for (bool muinv : {true, false})
    {
      auto A = muinv ? space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO)
                     : space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
      REQUIRE(A);
      A->Mult(x, y);
      ApplyReference(muinv, x, y_ref);
      INFO("muinv = " << muinv);
      CHECK_THAT(RelErr(), Matchers::WithinAbs(0.0, 1.0e-12));
    }
  }
  else
  {
    // K - ω² M + A2(ω), with the bulk terms in K and M, and the PML terms in A2(ω).
    auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto A2 = space_op.GetExtraSystemPMLMatrix(omega, Operator::DIAG_ZERO);
    REQUIRE(K);
    REQUIRE(M);
    REQUIRE(A2);
    K->Mult(x, y);
    M->Mult(x, t);
    y.Add(-omega * omega, t);
    A2->AddMult(x, y);
    ApplyReference(true, x, y_ref);
    ApplyReference(false, x, t);
    y_ref.Add(-omega * omega, t);
    CHECK_THAT(RelErr(), Matchers::WithinAbs(0.0, 1.0e-12));
  }
  BilinearForm::pa_order_threshold = pa_order_threshold;
  fem::DefaultIntegrationOrder::p_trial = p_trial;
}

TEST_CASE("SpaceOperator normal extruded PML matches the Cartesian PML",
          "[spaceoperator][pml][Serial][Parallel]")
{
  // The layer extruded along the normal +z of the top of a box, and the Cartesian PML block
  // of config["Domains"]["PML"] on the same elements, with the direction +z and the same
  // thickness and profile, give the same operators for an anisotropic, lossy background,
  // and so do the layer of the whole box, extruded per component, and the block with all of
  // its faces, edges, and corners, for the static stretch and the frequency-dependent one
  // with α > 0 at a complex frequency. After a uniform refinement, both measure the element
  // size of the layer on the mesh, half of t / N, for the same default σ_max.
  MPI_Comm comm = Mpi::World();
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  const bool box = GENERATE(false, true);
  const bool frequency_dependent = GENERATE(false, true);
  const bool refined = GENERATE(false, true);
  auto serial_mesh =
      std::make_unique<mfem::Mesh>(box ? MakeBox(type, 2, 2) : MakeBoxTop(type, 2, 2));
  json entry = {{"Attributes", {box ? 1 : 2}}};
  if (!refined)
  {
    entry["SigmaMax"] = 2.0;
  }
  auto extruded = MakePML(entry);
  ExtrudeEntry(serial_mesh, extruded);
  if (refined)
  {
    serial_mesh->UniformRefinement();
  }
  REQUIRE(extruded.extrusion.depth == (box ? config::PMLData::Extrusion::Depth::BOX
                                           : config::PMLData::Extrusion::Depth::NORMAL));
  config::PMLData cartesian;
  cartesian.attributes = extruded.attributes;
  cartesian.autodetect_geometry = box;
  cartesian.directions[5] = true;
  cartesian.thickness[5] = 0.5;
  if (!refined)
  {
    cartesian.sigma_auto = false;
    cartesian.sigma_max.fill(2.0);
  }
  for (auto *data : {&extruded, &cartesian})
  {
    data->reference_frequency = 1.5;
    data->kappa_max.fill(1.5);
    data->frequency_dependent = frequency_dependent;
    data->alpha_max.fill(frequency_dependent ? 0.3 : 0.0);
  }

  const int p_trial = fem::DefaultIntegrationOrder::p_trial;
  config::SolverData solver;
  solver.order = 2;
  solver.linear.mg_max_levels = 1;
  fem::DefaultIntegrationOrder::p_trial = solver.order;
  config::MaterialData material;
  material.attributes = {1, 2};
  material.mu_r.s = {1.0, 2.0, 1.5};
  material.epsilon_r.s = {2.0, 3.0, 4.0};
  material.tandelta.s = {1.0e-2, 2.0e-2, 3.0e-2};
  std::vector<std::unique_ptr<Mesh>> mesh[2];
  std::unique_ptr<SpaceOperator> space_op[2];
  for (int k = 0; k < 2; k++)
  {
    config::DomainData domains;
    domains.attributes = {1, 2};
    domains.materials = {material};
    domains.pml = {(k == 0) ? extruded : cartesian};
    mesh[k].push_back(std::make_unique<Mesh>(comm, *serial_mesh));
    space_op[k] =
        std::make_unique<SpaceOperator>(solver, domains, config::BoundaryData{},
                                        ProblemType::DRIVEN, Units(1.0, 1.0), mesh[k]);
    const auto &layer = space_op[k]->GetMaterialOp().GetPMLLayers().at(0);
    REQUIRE(layer.GetStretch().extruded == (k == 0 && !box));
    REQUIRE(layer.IsFrequencyDependent() == frequency_dependent);
    const auto &stretch = layer.GetStretch();
    for (int f = 0; f < 6; f++)
    {
      if (stretch.geometry.thickness[f] > 0.0)
      {
        INFO("k = " << k << ", f = " << f);
        CHECK(stretch.element_size[f] == Approx(refined ? 0.125 : 0.25));
      }
    }
  }
  // The stiffness and mass matrices, and the PML terms of A2(ω) of the frequency-dependent
  // stretch.
  for (int op = 0; op < (frequency_dependent ? 3 : 2); op++)
  {
    std::unique_ptr<ComplexOperator> A[2];
    for (int k = 0; k < 2; k++)
    {
      A[k] = (op == 0)
                 ? space_op[k]->GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO)
             : (op == 1)
                 ? space_op[k]->GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO)
                 : space_op[k]->GetExtraSystemPMLMatrix(1.3 + 0.05i, Operator::DIAG_ZERO);
      REQUIRE(A[k]);
    }
    ComplexVector x(A[0]->Width()), y(A[0]->Height()), y_ref(A[0]->Height());
    x.Real().Randomize(42 + Mpi::Rank(comm));
    x.Imag().Randomize(4242 + Mpi::Rank(comm));
    A[0]->Mult(x, y);
    A[1]->Mult(x, y_ref);
    const double norm = linalg::Norml2(comm, y_ref);
    y.Add(-1.0, y_ref);
    INFO("op = " << op);
    CHECK_THAT(linalg::Norml2(comm, y) / norm, Matchers::WithinAbs(0.0, 1.0e-12));
  }
  fem::DefaultIntegrationOrder::p_trial = p_trial;
}

namespace
{

// Driven configuration with config["Boundaries"]["PML"], of one entry or an array of them.
json MakeConfig(json pml, json boundaries = json::object())
{
  boundaries["PML"] = pml.is_array() ? pml : json::array({pml});
  return json{
      {"Problem", {{"Type", "Driven"}, {"Output", "test_output"}}},
      {"Model", {{"Mesh", "test.msh"}}},
      {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
      {"Boundaries", boundaries},
      {"Solver", {{"Driven", {{"MinFreq", 1.0}, {"MaxFreq", 3.0}, {"FreqStep", 1.0}}}}}};
}

}  // namespace

TEST_CASE("Config extruded PML", "[config][pml][Serial]")
{
  using Catch::Matchers::ContainsSubstring;
  const json pml = {{"Attributes", {4, 4}}, {"Thickness", 0.5}};

  // The extruded PML is the last PML block, with the default number of layers and core.
  {
    IoData iodata(MakeConfig(pml), false);
    REQUIRE(iodata.domains.pml.size() == 1);
    const auto &data = iodata.domains.pml[0];
    CHECK(data.IsExtruded());
    CHECK(data.attributes.empty());
    CHECK(data.extrusion.attributes == std::vector<int>{4});
    CHECK(data.extrusion.layers == 4);
    CHECK(!data.extrusion.has_center);
    CHECK(data.reference_frequency == 1.0);
    auto config = IoData::ConcretizeDefaults(iodata, MakeConfig(pml));
    CHECK(ValidateConfig(config).empty());
    CHECK(config["Boundaries"]["PML"][0]["Layers"] == 4);
    CHECK(config["Boundaries"]["PML"][0]["SigmaMax"] == "Auto");
    IoData iodata2(config, false);
    CHECK(iodata2.domains.pml.at(0).extrusion.layers == 4);
  }
  {
    auto config = MakeConfig(
        {{"Attributes", {4}}, {"Thickness", 0.5}, {"Layers", 3}, {"SigmaMax", 2.0}});
    CHECK(ValidateConfig(config).empty());
    IoData iodata(config, false);
    const auto &ex = iodata.domains.pml.at(0).extrusion;
    CHECK(ex.layers == 3);
    CHECK(iodata.domains.pml[0].sigma_max == std::array<double, 3>{2.0, 2.0, 2.0});
  }

  // Several entries are the last PML blocks, in order, each resolved with its own defaults,
  // and their surfaces must be disjoint.
  {
    const json entries = {{{"Attributes", {4}}, {"Thickness", 0.5}},
                          {{"Attributes", {5, 6}}, {"Thickness", 0.2}, {"Layers", 8}}};
    auto config = MakeConfig(entries);
    config["Domains"]["Materials"][0]["Attributes"] = {1, 2};
    config["Domains"]["PML"] = {{{"Attributes", {2}}}};
    CHECK(ValidateConfig(config).empty());
    IoData iodata(config, false);
    REQUIRE(iodata.domains.pml.size() == 3);
    CHECK(!iodata.domains.pml[0].IsExtruded());
    CHECK(iodata.domains.pml[1].extrusion.entry == 0);
    CHECK(iodata.domains.pml[2].extrusion.entry == 1);
    CHECK(iodata.domains.pml[2].extrusion.attributes == std::vector<int>{5, 6});
    const auto resolved = IoData::ConcretizeDefaults(iodata, config);
    CHECK(resolved["Boundaries"]["PML"][0]["Layers"] == 4);
    CHECK(resolved["Boundaries"]["PML"][1]["Layers"] == 8);
    CHECK(ValidateConfig(resolved).empty());
    auto overlap = entries;
    overlap[1]["Attributes"] = {4, 6};
    CHECK_THROWS_WITH(IoData(MakeConfig(overlap), false),
                      ContainsSubstring("is in more than one entry"));
  }

  // The entries form an array, and the extrusion is "Normal" or "Radial", with "Center"
  // only for the radial extrusion. The schema accepts "Center" with "Normal", so that the
  // error of the parser, which says why, reaches the user.
  {
    auto config = MakeConfig(pml);
    config["Boundaries"]["PML"] = pml;
    CHECK_FALSE(ValidateConfig(config).empty());
    CHECK_THROWS_WITH(IoData(config, false), ContainsSubstring("array of entries"));
  }
  for (const auto &[extrusion, center, valid] :
       std::vector<std::tuple<std::string, bool, bool>>{{"Normal", false, true},
                                                        {"Radial", false, true},
                                                        {"Radial", true, true},
                                                        {"Normal", true, false},
                                                        {"Ray", false, false}})
  {
    auto entry = pml;
    entry["Extrusion"] = extrusion;
    if (center)
    {
      entry["Center"] = {0.0, 0.0, 0.0};
    }
    INFO("Extrusion: " << extrusion << ", Center: " << center);
    CHECK(ValidateConfig(MakeConfig(entry)).empty() == (valid || center));
    if (valid)
    {
      IoData iodata(MakeConfig(entry), false);
      CHECK(iodata.domains.pml.at(0).extrusion.type ==
            ((extrusion == "Normal") ? config::PMLData::Extrusion::Type::NORMAL
                                     : config::PMLData::Extrusion::Type::RADIAL));
      const auto resolved = IoData::ConcretizeDefaults(iodata, MakeConfig(entry));
      CHECK(resolved["Boundaries"]["PML"][0]["Extrusion"] == extrusion);
      CHECK(ValidateConfig(resolved).empty());
    }
    else
    {
      CHECK_THROWS_WITH(IoData(MakeConfig(entry), false),
                        ContainsSubstring(center ? "is only used for the radial extrusion"
                                                 : "PML.Extrusion"));
    }
  }

  // The surface cannot have another boundary condition, the layer is refined uniformly or
  // by nonconformal adaptive refinement, and the stretch parameters are scalars.
  CHECK_THROWS_WITH(IoData(MakeConfig(pml, {{"Absorbing", {{"Attributes", {4}}}}}), false),
                    ContainsSubstring("cannot have another boundary condition"));
  {
    auto config = MakeConfig(pml);
    config["Model"]["Refinement"] = {{"UniformLevels", 1}};
    CHECK_NOTHROW(IoData(config, false));
    config["Model"]["Refinement"] = {{"MaxIts", 2}};
    CHECK_NOTHROW(IoData(config, false));
    config["Model"]["Refinement"]["Nonconformal"] = false;
    CHECK_THROWS_WITH(IoData(config, false),
                      ContainsSubstring("require nonconformal adaptive mesh refinement"));
    config["Problem"]["Type"] = "Transient";
    config["Model"].erase("Refinement");
    CHECK_THROWS(IoData(config, false));
  }
  auto with_array = pml;
  with_array["KappaMax"] = {1.0, 2.0, 1.0};
  CHECK_FALSE(ValidateConfig(MakeConfig(with_array)).empty());
  CHECK_THROWS_WITH(IoData(MakeConfig(with_array), false), ContainsSubstring("scalars"));
  CHECK_FALSE(ValidateConfig(MakeConfig({{"Attributes", {4}}})).empty());
}

TEST_CASE("Config extruded PML center", "[config][schema][pml][Serial]")
{
  // "Center" is a point or the two different corners of an axis-aligned box: the schema
  // enforces the form, and the parser gives the core as the box of the lower and upper
  // corners (a point if they are equal).
  using Catch::Matchers::ContainsSubstring;
  auto Config = [](json center)
  { return MakeConfig({{"Attributes", {4}}, {"Thickness", 0.5}, {"Center", center}}); };
  const json point = {1.0, 2.0, -3.0}, segment = {{0.0, 0.0, -1.0}, {0.0, 0.5, 1.0}};
  {
    // The corners in any order.
    IoData iodata(Config({{1.0, 0.0, 2.0}, {0.0, 1.0, -1.0}}), false);
    CHECK(iodata.domains.pml.at(0).extrusion.center ==
          std::array<std::array<double, 3>, 2>{{{0.0, 0.0, -1.0}, {1.0, 1.0, 2.0}}});
  }
  CHECK(ValidateConfig(Config(point)).empty());
  CHECK(ValidateConfig(Config(segment)).empty());
  CHECK(ValidateConfig(Config({0, 0, 1})).empty());
  for (const json &center :
       {json{0.0, 0.0}, json{0.0, 0.0, 0.0, 0.0}, json{{0.0, 0.0, -1.0}},
        json{{0.0, 0.0}, {0.0, 0.0, 1.0}}, json{{0.0, 0.0, -1.0}, {0.0, 0.0, 1.0, 2.0}},
        json{{0.0, 0.0, -1.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}},
        json{{0.0, 0.0, -1.0}, 1.0}, json{1.0, {0.0, 0.0, 1.0}},
        json{0.0, 0.0, {0.0, 0.0, 1.0}}, json("0 0 1"), json{"0", "0", "1"},
        json{{"0", "0", "-1"}, {"0", "0", "1"}},
        json{{{0.0, 0.0, -1.0}}, {{0.0, 0.0, 1.0}}},
        json{{{0.0, 0.0, -1.0}, {0.0, 0.0, 1.0}}}, json(1.0), json::object()})
  {
    INFO("Center: " << center.dump());
    CHECK_FALSE(ValidateConfig(Config(center)).empty());
    CHECK_THROWS_WITH(IoData(Config(center), false), ContainsSubstring("PML.Center"));
  }

  // The core of each form, and the resolved configuration keeps the form (or the default).
  for (const json &center : {point, segment})
  {
    IoData iodata(Config(center), false);
    const auto &ex = iodata.domains.pml.at(0).extrusion;
    CHECK(ex.has_center);
    CHECK(json(ex.center) == ((center == point) ? json{point, point} : segment));
    const auto resolved = IoData::ConcretizeDefaults(iodata, Config(center));
    CHECK(resolved["Boundaries"]["PML"][0]["Center"] == center);
    CHECK(ValidateConfig(resolved).empty());
    IoData iodata2(resolved, false);
    CHECK(iodata2.domains.pml.at(0).extrusion.center == ex.center);
  }
  {
    const auto config = MakeConfig({{"Attributes", {4}}, {"Thickness", 0.5}});
    IoData iodata(config, false);
    CHECK(!iodata.domains.pml.at(0).extrusion.has_center);
    CHECK(!IoData::ConcretizeDefaults(iodata, config)["Boundaries"]["PML"][0].contains(
        "Center"));
  }

  // A box with two equal corners is the point form.
  const json equal = {{1.0, 2.0, -3.0}, {1.0, 2.0, -3.0}};
  CHECK(ValidateConfig(Config(equal)).empty());
  CHECK_THROWS_WITH(IoData(Config(equal), false), ContainsSubstring("two equal corners"));
}

}  // namespace palace
