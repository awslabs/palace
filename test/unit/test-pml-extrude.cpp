// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <memory>
#include <string>
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
#include "linalg/operator.hpp"
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
// sphere (x ↦ x |x|_∞ / |x|).
mfem::Mesh MakeBall(mfem::Element::Type type, int n, int order)
{
  auto mesh = MakeBox(type, n, order);
  mesh.Transform(
      [](const mfem::Vector &x, mfem::Vector &y)
      {
        y = x;
        if (x.Norml2() > 0.0)
        {
          y *= x.Normlinf() / x.Norml2();
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
  // order, is at the radius 1 + ξ for its depth ξ = (k + z) t / N into the layer.
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  const int order = GENERATE(1, 2, 3);
  auto mesh = std::make_unique<mfem::Mesh>(MakeBall(type, 2, order));
  const int ne = mesh->GetNE(), nbe = mesh->GetNBE();
  auto pml = MakePML();
  mesh::ExtrudePML(mesh, pml, {});

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
  // those of the back of the layer, whose nodes are on the rays from the centre.
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  const int order = GENERATE(1, 2);
  auto mesh = std::make_unique<mfem::Mesh>(MakeBox(type, 2, order));
  const int ne = mesh->GetNE(), nbe = mesh->GetNBE();
  auto pml = MakePML();
  mesh::ExtrudePML(mesh, pml, {});
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
      // x = x⁰ (1 + t / |x⁰|) for the point x⁰ on the box.
      auto x = NodePosition(*mesh, be, j, true);
      x *= (x.Norml2() - 0.5) / x.Norml2();
      CHECK_THAT(x.Normlinf(), Matchers::WithinAbs(1.0, 1.0e-12));
    }
  }
}

TEST_CASE("ExtrudePML with a segment core", "[pml][geodata][Serial]")
{
  // With a segment core through the interface z = -0.5 between two materials of the box,
  // the rays of the faces beside the segment are normal to it, and the interface stays
  // planar through
  // the layer.
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  auto mesh = std::make_unique<mfem::Mesh>(MakeBox(type, 4, 2, true));
  const int ne = mesh->GetNE();
  auto pml = MakePML({{"Center", {{0.0, 0.0, -0.75}, {0.0, 0.0, 0.75}}}});
  mesh::ExtrudePML(mesh, pml, {});
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

  // With the default point core, the rays from the centre tilt the interface into a cone.
  mesh = std::make_unique<mfem::Mesh>(MakeBox(type, 4, 2, true));
  pml = MakePML();
  mesh::ExtrudePML(mesh, pml, {});
  z_above = mfem::infinity();
  for (int e = ne; e < mesh->GetNE(); e++)
  {
    for (int j = 0; j < mesh->GetNodes()->FESpace()->GetFE(e)->GetDof(); j++)
    {
      if (mesh->GetAttribute(e) == 4)
      {
        z_above = std::min(z_above, NodePosition(*mesh, e, j)(2));
      }
    }
  }
  CHECK(z_above < -0.5 - 1.0e-2);
}

TEST_CASE("ExtrudePML rejects invalid surfaces", "[pml][geodata][Serial]")
{
  using Catch::Matchers::ContainsSubstring;
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  auto Extrude =
      [&](json extra, const config::BoundaryData &boundaries = {}, bool interior = false)
  {
    auto mesh = std::make_unique<mfem::Mesh>(MakeBox(type, 2, 1, true));
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
    mesh::ExtrudePML(mesh, pml, boundaries);
  };

  // A core outside of the box: the face x = 1 faces away from it.
  CHECK_THROWS_WITH(Extrude({{"Center", {3.0, 0.0, 0.0}}}),
                    ContainsSubstring("does not face its closest point on the point"));
  // A core on the surface.
  CHECK_THROWS_WITH(Extrude({{"Center", {1.0, 0.0, 0.0}}}),
                    ContainsSubstring("does not face its closest point on the point"));
  // An interior boundary.
  CHECK_THROWS_WITH(Extrude(json::object(), {}, true),
                    ContainsSubstring("is not on the exterior boundary"));
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
    mesh::ExtrudePML(mesh, pml, boundaries);
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
  mesh::ExtrudePML(mesh, pml, {});
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

TEST_CASE("Extruded PML depth", "[pml][geodata][Serial][Parallel]")
{
  // The depth of the nodes of the layer, identified by their rays from the surface on all
  // processes, is the depth of the extrusion: on the sphere, the radius less one.
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  auto smesh = std::make_unique<mfem::Mesh>(MakeBall(type, 2, 2));
  const int ne = smesh->GetNE();
  auto pml = MakePML();
  mesh::ExtrudePML(smesh, pml, {});
  // The outer layer is on the last process, without the surface.
  const int size = Mpi::Size(Mpi::World());
  std::vector<int> partitioning(smesh->GetNE());
  for (int e = 0; e < smesh->GetNE(); e++)
  {
    partitioning[e] = (e >= ne && (e - ne) % 2 == 1) ? size - 1 : e % std::max(size - 1, 1);
  }
  mfem::ParMesh pmesh(Mpi::World(), *smesh, partitioning.data());
  const auto depth = mesh::GetExtrudedPMLDepth(pmesh, pml);
  double err = 0.0, depth_max = 0.0;
  mfem::Array<int> dofs;
  for (int e = 0; e < pmesh.GetNE(); e++)
  {
    const bool in_layer = (pmesh.GetAttribute(e) == 2);
    depth->FESpace()->GetElementDofs(e, dofs);
    for (int j = 0; j < dofs.Size(); j++)
    {
      const double r = NodePosition(pmesh, e, j).Norml2();
      err = std::max(err, std::abs((*depth)(dofs[j]) - (in_layer ? r - 1.0 : 0.0)));
      depth_max = std::max(depth_max, (*depth)(dofs[j]));
    }
  }
  Mpi::GlobalMax(1, &err, pmesh.GetComm());
  Mpi::GlobalMax(1, &depth_max, pmesh.GetComm());
  CHECK(err < 1.0e-12);
  CHECK_THAT(depth_max, Matchers::WithinAbs(0.5, 1.0e-12));

  // On a refined box, the nodes of the layer are not on the rays of the surface.
  smesh = std::make_unique<mfem::Mesh>(MakeBox(type, 2, 1));
  pml = MakePML();
  mesh::ExtrudePML(smesh, pml, {});
  mfem::ParMesh pmesh_ref(Mpi::World(), *smesh);
  pmesh_ref.UniformRefinement();
  CHECK_THROWS_WITH(mesh::GetExtrudedPMLDepth(pmesh_ref, pml),
                    Matchers::ContainsSubstring("cannot be refined"));
}

namespace
{

// Extruded stretch around a point core at the origin, with the profile of a layer of
// thickness 0.5.
pml::Stretch MakeExtrudedStretch(bool frequency_dependent)
{
  pml::Stretch p;
  p.extruded = true;
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
  auto s = [&](double x)
  {
    const double shape = std::pow(x / p.geometry.thickness[0], p.order);
    return 1.0 + (p.kappa_max[0] - 1.0) * shape +
           p.sigma_max[0] * shape / (p.alpha_max[0] * shape + 1i * omega);
  };
  constexpr int n = 2000;
  std::complex<double> integral = 0.0;
  for (int i = 0; i < n; i++)
  {
    const double a = xi * i / n, b = xi * (i + 1) / n;
    integral += (b - a) / 6.0 * (s(a) + 4.0 * s(0.5 * (a + b)) + s(b));
  }
  return {s(xi), integral};
}

// Requested part of c μ̃⁻¹ or c ε̃ from the device helpers, for an isotropic or anisotropic
// background.
std::array<double, 9> EvalCoeffAt(const pml::Stretch &p, const pml::Background &b,
                                  pml::TensorPart part, std::complex<double> c,
                                  std::complex<double> omega, bool muinv,
                                  const std::array<double, 3> &x, double xi,
                                  const std::array<double, 3> &grad_xi)
{
  pml::ContextHeader header;
  header.part = part;
  (muinv ? header.c_muinv : header.c_eps) = c;
  header.omega = omega;
  const auto ctx = pml::PackContext(header, p, {0}, {b});
  std::array<double, 9> coeff;
  REQUIRE(PMLCoeffAt(ctx.data(), 1, x.data(), xi, grad_xi.data(),
                     muinv ? coeff.data() : nullptr, muinv ? nullptr : coeff.data()));
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
    case pml::TensorPart::ABS:
      return c.real() * std::abs(t);
    case pml::TensorPart::STRETCH_MAGNITUDE_REAL:
      return (c * t_mag).real();
    case pml::TensorPart::STRETCH_MAGNITUDE_IMAG:
      return (c * t_mag).imag();
  }
  return 0.0;
}

constexpr std::array<pml::TensorPart, 5> all_parts = {
    pml::TensorPart::REAL, pml::TensorPart::IMAG, pml::TensorPart::ABS,
    pml::TensorPart::STRETCH_MAGNITUDE_REAL, pml::TensorPart::STRETCH_MAGNITUDE_IMAG};

}  // namespace

TEST_CASE("Extruded PML stretch", "[pml][Serial]")
{
  const bool frequency_dependent = GENERATE(false, true);
  const auto p = MakeExtrudedStretch(frequency_dependent);
  const std::complex<double> omega = 1.3 + 0.05i, c = 0.7 - 0.2i;

  SECTION("Identity at the interface")
  {
    const std::array<double, 3> x = {0.3, -0.5, 1.1}, grad_xi = {0.2, 0.7, -0.1};
    for (bool magnitude : {false, true})
    {
      const auto J = pml::EvaluateJacobian(p, x, 0.0, grad_xi, omega, magnitude);
      for (int k = 0; k < 9; k++)
      {
        CHECK(J[k] == std::complex<double>((k % 4 == 0) ? 1.0 : 0.0, 0.0));
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

      // The principal stretch factors along and across the rays, from which the smoother
      // masks of the multigrid preconditioner are selected.
      const auto f = pml::EvaluateExtrudedStretch(p, x, xi, omega);
      CHECK(std::abs(f[0] - s) <= 1.0e-9 * std::abs(s));
      CHECK(std::abs(f[1] - q) <= 1.0e-9 * std::abs(q));
      for (bool muinv : {true, false})
      {
        for (auto part : all_parts)
        {
          const auto coeff = EvalCoeffAt(p, b, part, c, omega, muinv, x, xi, dir);
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

  SECTION("Cartesian PML on a face of a box")
  {
    // For rays normal to a face (d ∇ξᵀ = n nᵀ, ∇d = 0), J̃ = I + (s - 1) n nᵀ is the
    // diagonal stretch S of the Cartesian PML of the face, and the tensors are the same for
    // an anisotropic, lossy background.
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
      for (bool muinv : {true, false})
      {
        for (auto part : all_parts)
        {
          pml::ContextHeader header;
          header.part = part;
          (muinv ? header.c_muinv : header.c_eps) = c;
          header.omega = omega;
          const auto ctx_s = pml::PackContext(header, cartesian, {0}, {b});
          const auto ctx_j = pml::PackContext(header, p, {0}, {b});
          std::array<double, 9> coeff_s, coeff_j;
          muinv ? PMLMuInvCoeff(ctx_s.data(), 1, x.data(), coeff_s.data())
                : PMLEpsCoeff(ctx_s.data(), 1, x.data(), coeff_s.data());
          const bool magnitude = (part == pml::TensorPart::STRETCH_MAGNITUDE_REAL ||
                                  part == pml::TensorPart::STRETCH_MAGNITUDE_IMAG);
          const auto s_part = magnitude ? std::complex<double>(std::abs(s)) : s;
          CeedScalar J_re[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, s_part.real()},
                     J_im[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, s_part.imag()};
          const CeedIntScalar *bg = PMLBackgroundData(ctx_j.data(), 1);
          muinv ? PMLMuInvCoeffJacobian(ctx_j.data(), bg, J_re, J_im, coeff_j.data())
                : PMLEpsCoeffJacobian(ctx_j.data(), bg, J_re, J_im, coeff_j.data());
          for (int k = 0; k < 9; k++)
          {
            CHECK_THAT(coeff_j[k], Matchers::WithinAbs(coeff_s[k], 1.0e-12));
          }
        }
      }
    }
  }
}

namespace
{

// The tensors of the extruded stretch at a quadrature point, evaluated independently of the
// device helpers from the depth of the mesh nodes, for a point core at the origin and an
// isotropic background.
class ReferenceExtrudedTensor : public mfem::MatrixCoefficient
{
private:
  const pml::Stretch &p;
  const mfem::GridFunction &depth;
  double mu_inv;
  std::complex<double> eps;
  bool muinv, imag;

public:
  ReferenceExtrudedTensor(const pml::Stretch &p, const mfem::GridFunction &depth,
                          double mu_inv, std::complex<double> eps, bool muinv, bool imag)
    : mfem::MatrixCoefficient(3), p(p), depth(depth), mu_inv(mu_inv), eps(eps),
      muinv(muinv), imag(imag)
  {
  }

  void Eval(mfem::DenseMatrix &K, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    T.SetIntPoint(&ip);
    mfem::Vector x(3), grad_xi(3);
    T.Transform(ip, x);
    const double xi = depth.GetValue(T, ip);
    depth.GetGradient(T, grad_xi);
    const auto [s, xi_t] = ReferenceProfile(p, xi, 0.0);
    const double r = x.Norml2();
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
  // Ball of hexahedra with an extruded spherical layer: the libCEED operators of the PML
  // terms, from the depth and its gradient in the geometry factor data, match MFEM
  // integrators with the tensors evaluated from the depth of the mesh nodes.
  MPI_Comm comm = Mpi::World();
  auto serial_mesh =
      std::make_unique<mfem::Mesh>(MakeBall(mfem::Element::HEXAHEDRON, 2, 2));
  config::DomainData domains;
  domains.pml = {MakePML()};
  domains.pml[0].reference_frequency = 1.5;
  domains.pml[0].kappa_max.fill(1.5);
  mesh::ExtrudePML(serial_mesh, domains.pml[0], {});
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(comm, *serial_mesh));

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
  const auto *depth = mesh[0]->GetPMLDepth();
  REQUIRE(depth);

  const auto &pfes = space_op.GetNDSpace().Get();
  const auto &ir = mfem::IntRules.Get(
      mfem::Geometry::CUBE,
      fem::DefaultIntegrationOrder::Get(mesh[0]->Get(), mfem::Geometry::CUBE));
  const std::complex<double> eps = 2.5;
  auto AssembleReference = [&](bool muinv, bool imag)
  {
    mfem::ParBilinearForm a(const_cast<mfem::ParFiniteElementSpace *>(&pfes));
    mfem::Array<int> pml_marker({0, 1}), bulk_marker({1, 0});
    ReferenceExtrudedTensor pml_coef(p, *depth, 1.0 / 1.5, eps, muinv, imag);
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

  for (bool muinv : {true, false})
  {
    auto Ar = AssembleReference(muinv, false), Ai = AssembleReference(muinv, true);
    auto A = muinv ? space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO)
                   : space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    REQUIRE(A);
    ComplexVector x(A->Width()), y(A->Height()), y_ref(A->Height());
    x.Real().Randomize(42 + Mpi::Rank(comm));
    x.Imag().Randomize(4242 + Mpi::Rank(comm));
    A->Mult(x, y);
    Vector t(A->Height());
    Ar->Mult(x.Real(), y_ref.Real());
    Ai->Mult(x.Imag(), t);
    y_ref.Real().Add(-1.0, t);
    Ar->Mult(x.Imag(), y_ref.Imag());
    Ai->Mult(x.Real(), t);
    y_ref.Imag().Add(1.0, t);
    const double norm = linalg::Norml2(comm, y_ref);
    y.Add(-1.0, y_ref);
    INFO("muinv = " << muinv);
    CHECK_THAT(linalg::Norml2(comm, y) / norm, Matchers::WithinAbs(0.0, 1.0e-12));
  }
  BilinearForm::pa_order_threshold = pa_order_threshold;
  fem::DefaultIntegrationOrder::p_trial = p_trial;
}

TEST_CASE("SpaceOperator smoother masks of an extruded PML",
          "[spaceoperator][pml][Serial][Parallel]")
{
  // Ball of hexahedra with an extruded spherical layer: the multigrid smoothers do not
  // update the unknowns of the strongly stretched extruded elements, selected with the
  // stretch factors along and across the rays at the depth of their vertices (the stretch
  // of the Cartesian PML is one everywhere for an extruded layer). The masks for
  // increasing thresholds are nested, and contained in the PML subdomain.
  MPI_Comm comm = Mpi::World();
  const int pa_order_threshold = BilinearForm::pa_order_threshold;
  const int p_trial = fem::DefaultIntegrationOrder::p_trial;
  config::SolverData solver;
  solver.order = 2;
  solver.linear.mg_max_levels = 2;
  solver.linear.mg_coarsening = MultigridCoarsening::LINEAR;
  solver.linear.pml_subdomain_solver = PMLSubdomainSolver::DIRECT;
  solver.linear.pc_mat_real = false;
  solver.linear.pc_mat_shifted = 0;
  fem::DefaultIntegrationOrder::p_trial = solver.order;
  config::MaterialData material;
  material.attributes = {1, 2};
  config::BoundaryData boundaries;
  Units units(1.0, 1.0);

  constexpr double omega = 1.5;
  auto Select = [&](double smoother_threshold)
  {
    // The depth of the nodes into the layer is set on the mesh by the SpaceOperator, so
    // each one has its own mesh.
    auto serial_mesh =
        std::make_unique<mfem::Mesh>(MakeBall(mfem::Element::HEXAHEDRON, 2, 2));
    config::DomainData domains;
    domains.pml = {MakePML()};
    domains.pml[0].reference_frequency = 1.5;
    mesh::ExtrudePML(serial_mesh, domains.pml[0], {});
    domains.attributes = {1, 2};
    domains.materials = {material};
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(comm, *serial_mesh));
    solver.linear.pml_smoother_threshold = smoother_threshold;
    SpaceOperator space_op(solver, domains, boundaries, ProblemType::DRIVEN, units, mesh);
    REQUIRE(space_op.GetMaterialOp().GetPMLLayers().at(0).GetStretch().extruded);
    auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(
        1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i, omega);
    const auto *mg = dynamic_cast<const ComplexMultigridOperator *>(P.get());
    REQUIRE(mg);
    REQUIRE(mg->GetNumLevels() == 2);
    return std::array{mg->GetPMLTrueDofs(), mg->GetSmootherMask(1),
                      mg->GetAuxiliarySmootherMask(1)};
  };
  auto GlobalSize = [&](const std::vector<int> &v)
  {
    long long n = static_cast<long long>(v.size());
    Mpi::GlobalSum(1, &n, comm);
    return n;
  };
  const auto all = Select(1.0), stretched = Select(1.01), none = Select(1.0e6);
  CHECK(GlobalSize(all[0]) > 0);
  CHECK(all[1] == all[0]);
  CHECK(stretched[0] == all[0]);
  CHECK(std::ranges::includes(all[1], stretched[1]));
  CHECK(GlobalSize(stretched[1]) > 0);
  CHECK(GlobalSize(stretched[2]) > 0);
  CHECK(GlobalSize(none[1]) == 0);
  CHECK(GlobalSize(none[2]) == 0);
  BilinearForm::pa_order_threshold = pa_order_threshold;
  fem::DefaultIntegrationOrder::p_trial = p_trial;
}

namespace
{

// Driven configuration with config["Boundaries"]["PML"].
json MakeConfig(json pml, json boundaries = json::object())
{
  boundaries["PML"] = pml;
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
    CHECK(config["Boundaries"]["PML"]["Layers"] == 4);
    CHECK(config["Boundaries"]["PML"]["SigmaMax"] == "Auto");
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

  // The surface cannot have another boundary condition, the layer cannot be refined, and
  // the stretch parameters are scalars.
  CHECK_THROWS_WITH(IoData(MakeConfig(pml, {{"Absorbing", {{"Attributes", {4}}}}}), false),
                    ContainsSubstring("cannot have another boundary condition"));
  {
    auto config = MakeConfig(pml);
    config["Model"]["Refinement"] = {{"UniformLevels", 1}};
    CHECK_THROWS_WITH(IoData(config, false), ContainsSubstring("refinement"));
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
  // "Center" is a point or a segment between two different points: the schema enforces the
  // form, and the parser gives the core as the segment between the two points (a point if
  // they are equal).
  using Catch::Matchers::ContainsSubstring;
  auto Config = [](json center)
  { return MakeConfig({{"Attributes", {4}}, {"Thickness", 0.5}, {"Center", center}}); };
  const json point = {1.0, 2.0, -3.0}, segment = {{0.0, 0.0, -1.0}, {0.0, 0.5, 1.0}};
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
    CHECK(resolved["Boundaries"]["PML"]["Center"] == center);
    CHECK(ValidateConfig(resolved).empty());
    IoData iodata2(resolved, false);
    CHECK(iodata2.domains.pml.at(0).extrusion.center == ex.center);
  }
  {
    const auto config = MakeConfig({{"Attributes", {4}}, {"Thickness", 0.5}});
    IoData iodata(config, false);
    CHECK(!iodata.domains.pml.at(0).extrusion.has_center);
    CHECK(!IoData::ConcretizeDefaults(iodata, config)["Boundaries"]["PML"].contains(
        "Center"));
  }

  // A segment between two equal points is the point form.
  const json equal = {{1.0, 2.0, -3.0}, {1.0, 2.0, -3.0}};
  CHECK(ValidateConfig(Config(equal)).empty());
  CHECK_THROWS_WITH(IoData(Config(equal), false), ContainsSubstring("two equal points"));
}

}  // namespace palace
