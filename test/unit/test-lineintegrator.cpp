// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <functional>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "fem/lineintegrator.hpp"
#include "utils/communication.hpp"

namespace palace
{
using namespace Catch::Matchers;

namespace
{

using VecFunc = std::function<void(const mfem::Vector &, mfem::Vector &)>;

// Reference value of ∫ F · dl along a polyline, using a high-order Gauss rule on each
// segment (exact for the polynomial fields used here).
double ReferenceLineIntegral(const VecFunc &F, const std::vector<mfem::Vector> &path,
                             bool closed, int vdim)
{
  const auto &ir = mfem::IntRules.Get(mfem::Geometry::SEGMENT, 20);
  const int sdim = path[0].Size();
  const int n = static_cast<int>(path.size()) - (closed ? 0 : 1);
  double result = 0.0;
  mfem::Vector x(sdim), f(vdim);
  for (int s = 0; s < n; s++)
  {
    const auto &p0 = path[s];
    const auto &p1 = path[(s + 1) % path.size()];
    for (int i = 0; i < ir.GetNPoints(); i++)
    {
      const auto &ip = ir.IntPoint(i);
      for (int d = 0; d < sdim; d++)
      {
        x(d) = p0(d) + ip.x * (p1(d) - p0(d));
      }
      F(x, f);
      for (int d = 0; d < std::min(sdim, vdim); d++)
      {
        result += ip.weight * f(d) * (p1(d) - p0(d));
      }
    }
  }
  return result;
}

mfem::Vector Point(std::initializer_list<double> x)
{
  mfem::Vector v(static_cast<int>(x.size()));
  int i = 0;
  for (auto xi : x)
  {
    v(i++) = xi;
  }
  return v;
}

// Locally refine the elements whose centers satisfy the predicate, to mimic adaptive
// refinement towards a point (conforming bisection or nonconforming refinement).
void RefineNear(mfem::ParMesh &mesh, const mfem::Vector &c, double r, int levels,
                bool nonconforming)
{
  mfem::Vector x(mesh.SpaceDimension());
  for (int l = 0; l < levels; l++)
  {
    mfem::Array<int> marked;
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      auto &T = *mesh.GetElementTransformation(e);
      T.Transform(mfem::Geometries.GetCenter(T.GetGeometryType()), x);
      if (x.DistanceTo(c) < r)
      {
        marked.Append(e);
      }
    }
    mesh.GeneralRefinement(marked, nonconforming ? 1 : 0);
    r *= 0.6;
  }
}

}  // namespace

TEST_CASE("LineIntegrator 2D exact for discrete fields",
          "[lineintegrator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const bool nonconforming = GENERATE(false, true);
  const auto type = GENERATE(mfem::Element::TRIANGLE, mfem::Element::QUADRILATERAL);
  if (type == mfem::Element::QUADRILATERAL && !nonconforming)
  {
    return;  // Local conforming refinement is only available for simplices
  }
  CAPTURE(nonconforming, type);

  auto smesh = mfem::Mesh::MakeCartesian2D(6, 4, type, false, 3.0, 2.0);
  if (nonconforming)
  {
    smesh.EnsureNCMesh(true);
  }
  mfem::ParMesh mesh(comm, smesh);

  // Graded refinement towards a point on the path, as AMR does near a metal edge.
  RefineNear(mesh, Point({1.0, 1.0}), 0.6, 6, nonconforming);

  // Field components are quadratic polynomials, contained in the order 3 Nédélec space.
  const int order = 3;
  mfem::ND_FECollection fec(order, 2);
  mfem::ParFiniteElementSpace fespace(&mesh, &fec);
  VecFunc F = [](const mfem::Vector &x, mfem::Vector &f)
  {
    f(0) = 1.0 + x(0) * x(0) + 0.5 * x(1);
    f(1) = x(0) * x(1) - x(1) * x(1) + 2.0;
  };
  mfem::VectorFunctionCoefficient coeff(2, F);
  mfem::ParGridFunction U(&fespace);
  U.ProjectCoefficient(coeff);

  auto Check = [&](const std::vector<mfem::Vector> &path, bool closed)
  {
    fem::LineIntegrator integrator(mesh, path, closed, 2 * order);
    const double ref = ReferenceLineIntegral(F, path, closed, 2);
    CAPTURE(ref);
    CHECK_THAT(integrator.Integrate(U), WithinAbs(ref, 1.0e-10 * (1.0 + std::abs(ref))));
    CHECK_THAT(integrator.GetUncoveredLength(), WithinAbs(0.0, 1.0e-12));
  };

  SECTION("Oblique path through the refined region")
  {
    Check({Point({0.13, 0.27}), Point({1.0, 1.0}), Point({2.71, 1.83})}, false);
  }

  SECTION("Path along element edges")
  {
    // y = 1 and x = 1.5 are grid lines of the coarse mesh, so every piece of these
    // segments lies on an edge shared by two elements: it must be counted exactly once.
    Check({Point({0.0, 1.0}), Point({1.5, 1.0}), Point({1.5, 2.0})}, false);
  }

  SECTION("Closed path")
  {
    Check({Point({0.3, 0.3}), Point({2.4, 0.6}), Point({1.1, 1.7})}, true);
  }

  SECTION("Path leaving the mesh")
  {
    // Only the portion inside the mesh contributes.
    std::vector<mfem::Vector> path = {Point({-1.0, 0.5}), Point({1.0, 0.5})};
    fem::LineIntegrator integrator(mesh, path, false, 2 * order);
    const double ref =
        ReferenceLineIntegral(F, {Point({0.0, 0.5}), Point({1.0, 0.5})}, false, 2);
    CHECK_THAT(integrator.Integrate(U), WithinAbs(ref, 1.0e-10 * (1.0 + std::abs(ref))));
    CHECK_THAT(integrator.GetUncoveredLength(), WithinAbs(1.0, 1.0e-10));
    CHECK_THAT(integrator.GetPathLength(), WithinAbs(2.0, 1.0e-14));
  }
}

TEST_CASE("LineIntegrator on a 3D boundary submesh", "[lineintegrator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  CAPTURE(type);

  // MakeCartesian3D boundary attribute 5 is the x = 0 face.
  auto smesh = mfem::Mesh::MakeCartesian3D(4, 4, 4, type, 1.0, 1.0, 1.0);
  smesh.EnsureNCMesh(true);
  mfem::ParMesh mesh(comm, smesh);
  RefineNear(mesh, Point({0.0, 0.5, 0.5}), 0.4, 3, true);
  mfem::Array<int> attr(1);
  attr[0] = 5;
  auto submesh = mfem::ParSubMesh::CreateFromBoundary(mesh, attr);

  // Linear tangential field, contained in the order 2 Nédélec space on the surface.
  const int order = 2;
  mfem::ND_FECollection fec(order, 2);
  mfem::ParFiniteElementSpace fespace(&submesh, &fec);
  VecFunc F = [](const mfem::Vector &x, mfem::Vector &f)
  {
    f(0) = 7.0;  // Normal component, ignored
    f(1) = 1.0 + x(1) + 2.0 * x(2);
    f(2) = x(1) - 2.0 * x(2) + 0.5;
  };
  mfem::VectorFunctionCoefficient coeff(3, F);
  mfem::ParGridFunction U(&fespace);
  U.ProjectCoefficient(coeff);

  for (const auto &path :
       {std::vector<mfem::Vector>{Point({0.0, 0.1, 0.2}), Point({0.0, 0.5, 0.5}),
                                  Point({0.0, 0.9, 0.7})},
        std::vector<mfem::Vector>{Point({0.0, 0.5, 0.0}), Point({0.0, 0.5, 1.0})}})
  {
    fem::LineIntegrator integrator(submesh, path, false, 2 * order);
    const double ref = ReferenceLineIntegral(F, path, false, 3);
    CHECK_THAT(integrator.Integrate(U), WithinAbs(ref, 1.0e-10 * (1.0 + std::abs(ref))));
    CHECK_THAT(integrator.GetUncoveredLength(), WithinAbs(0.0, 1.0e-12));
  }
}

TEST_CASE("LineIntegrator exact for discrete gradient fields",
          "[lineintegrator][Serial][Parallel]")
{
  // For E = ∇φ with φ a random H1 finite element function, E lies in the Nédélec space and
  // is discontinuous across element boundaries, while its line integral along any path is
  // exactly φ(p1) - φ(p0). A single global Gauss rule over the segment would not resolve
  // the kinks of the integrand, but the element-aware quadrature must be exact.
  MPI_Comm comm = Mpi::World();
  const bool nonconforming = GENERATE(false, true);
  CAPTURE(nonconforming);
  auto smesh = mfem::Mesh::MakeCartesian2D(5, 4, mfem::Element::TRIANGLE, false, 3.0, 2.0);
  if (nonconforming)
  {
    smesh.EnsureNCMesh(true);
  }
  mfem::ParMesh mesh(comm, smesh);
  RefineNear(mesh, Point({1.2, 1.0}), 0.6, 8, nonconforming);

  const int order = 2;
  mfem::H1_FECollection h1_fec(order, 2);
  mfem::ND_FECollection nd_fec(order, 2);
  mfem::ParFiniteElementSpace h1_fespace(&mesh, &h1_fec), nd_fespace(&mesh, &nd_fec);
  mfem::Vector phi_t(h1_fespace.GetTrueVSize()), E_t(nd_fespace.GetTrueVSize());
  phi_t.Randomize(1 + Mpi::Rank(comm));
  mfem::ParGridFunction phi(&h1_fespace), E(&nd_fespace);
  phi.SetFromTrueDofs(phi_t);
  {
    mfem::ParDiscreteLinearOperator grad(&h1_fespace, &nd_fespace);
    grad.AddDomainInterpolator(new mfem::GradientInterpolator);
    grad.Assemble();
    grad.Finalize();
    std::unique_ptr<mfem::HypreParMatrix> G(grad.ParallelAssemble());
    G->Mult(phi_t, E_t);
  }
  E.SetFromTrueDofs(E_t);

  // Evaluate φ at the given points (collective).
  auto Phi = [&](const std::vector<mfem::Vector> &points)
  {
    mfem::DenseMatrix pts(2, static_cast<int>(points.size()));
    for (int i = 0; i < pts.Width(); i++)
    {
      pts(0, i) = points[i](0);
      pts(1, i) = points[i](1);
    }
    mfem::Array<int> elems;
    mfem::Array<mfem::IntegrationPoint> ips;
    mesh.FindPoints(pts, elems, ips, false);
    std::vector<double> vals(points.size(), 0.0);
    for (int i = 0; i < pts.Width(); i++)
    {
      if (elems[i] >= 0)
      {
        vals[i] = phi.GetValue(elems[i], ips[i]);
      }
    }
    Mpi::GlobalSum(static_cast<int>(vals.size()), vals.data(), comm);
    return vals;
  };

  const std::vector<mfem::Vector> path = {Point({0.17, 0.31}), Point({1.2, 1.0}),
                                          Point({1.2, 1.5}), Point({2.83, 1.77})};
  const auto vals = Phi({path.front(), path.back()});
  const double ref = vals[1] - vals[0];
  {
    fem::LineIntegrator integrator(mesh, path, false, 2 * order);
    CAPTURE(ref, integrator.GetNumSubIntervals());
    CHECK_THAT(integrator.Integrate(E), WithinAbs(ref, 1.0e-10 * (1.0 + std::abs(ref))));
  }
  {
    fem::LineIntegrator integrator(mesh, path, true, 2 * order);
    CHECK_THAT(integrator.Integrate(E), WithinAbs(0.0, 1.0e-10 * (1.0 + std::abs(ref))));
  }
}

TEST_CASE("LineIntegrator on curved elements", "[lineintegrator][Serial][Parallel]")
{
  // Quarter annulus a < r < b with a curved (cubic) mesh. Straight paths ending next to
  // the curved boundaries cross curved elements and end between mesh vertices, where the
  // straight chords of the boundary elements would cut off part of the path. The line
  // integral of E = ∇φ must still be φ(p1) - φ(p0).
  MPI_Comm comm = Mpi::World();
  const auto type = GENERATE(mfem::Element::TRIANGLE, mfem::Element::QUADRILATERAL);
  CAPTURE(type);
  const double a = 0.5, b = 2.0;
  auto smesh = mfem::Mesh::MakeCartesian2D(3, 8, type, false, 1.0, 1.0);
  smesh.SetCurvature(3);
  smesh.Transform(
      [&](const mfem::Vector &u, mfem::Vector &x)
      {
        const double r = a + (b - a) * u(0), theta = 0.5 * M_PI * u(1);
        x(0) = r * std::cos(theta);
        x(1) = r * std::sin(theta);
      });
  mfem::ParMesh mesh(comm, smesh);

  const int order = 2;
  mfem::H1_FECollection h1_fec(order, 2);
  mfem::ND_FECollection nd_fec(order, 2);
  mfem::ParFiniteElementSpace h1_fespace(&mesh, &h1_fec), nd_fespace(&mesh, &nd_fec);
  mfem::FunctionCoefficient phi_coeff(
      [](const mfem::Vector &x) { return std::sin(1.3 * x(0)) * std::cosh(0.7 * x(1)); });
  mfem::ParGridFunction phi(&h1_fespace), E(&nd_fespace);
  phi.ProjectCoefficient(phi_coeff);
  {
    mfem::ParDiscreteLinearOperator grad(&h1_fespace, &nd_fespace);
    grad.AddDomainInterpolator(new mfem::GradientInterpolator);
    grad.Assemble();
    grad.Finalize();
    grad.Mult(phi, E);
  }

  // Evaluate φ at a point (collective), searching all elements with the exact curved
  // geometry (points on the curved boundary between vertices lie outside the straight
  // chords used by a bounding-box search).
  auto Phi = [&](const mfem::Vector &x)
  {
    mfem::InverseElementTransformation inv;
    inv.SetSolverType(mfem::InverseElementTransformation::NewtonElementProject);
    mfem::IntegrationPoint ip;
    mfem::Vector y(2);
    double val = 0.0;
    int found = Mpi::Size(comm);
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      auto &T = *mesh.GetElementTransformation(e);
      inv.SetTransformation(T);
      if (inv.Transform(x, ip) == mfem::InverseElementTransformation::Unknown)
      {
        continue;
      }
      T.Transform(ip, y);
      if (y.DistanceTo(x) < 1.0e-12)
      {
        val = phi.GetValue(e, ip);
        found = Mpi::Rank(comm);
        break;
      }
    }
    int owner = found;
    Mpi::GlobalMin(1, &owner, comm);
    REQUIRE(owner < Mpi::Size(comm));
    val = (owner == Mpi::Rank(comm)) ? val : 0.0;
    Mpi::GlobalSum(1, &val, comm);
    return val;
  };

  // Endpoints just inside the mesh, which only approximates the circles between nodes, but
  // well within the gap between the curved boundary and its straight chords.
  const double th0 = 0.37 * 0.5 * M_PI, th1 = 0.81 * 0.5 * M_PI;
  const double ra = a + 1.0e-3, rb = b - 1.0e-3;
  for (const auto &path :
       {std::vector<mfem::Vector>{Point({ra * std::cos(th0), ra * std::sin(th0)}),
                                  Point({rb * std::cos(th0), rb * std::sin(th0)})},
        std::vector<mfem::Vector>{Point({rb * std::cos(th0), rb * std::sin(th0)}),
                                  Point({rb * std::cos(th1), rb * std::sin(th1)})}})
  {
    fem::LineIntegrator integrator(mesh, path, false, 2 * (order + 3));
    const double ref = Phi(path.back()) - Phi(path.front());
    CAPTURE(ref);
    CHECK_THAT(integrator.Integrate(E), WithinAbs(ref, 1.0e-8 * (1.0 + std::abs(ref))));
    CHECK_THAT(integrator.GetUncoveredLength(), WithinAbs(0.0, 1.0e-12));
  }
}

}  // namespace palace
