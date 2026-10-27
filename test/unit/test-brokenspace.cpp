// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "fem/bilinearform.hpp"
#include "fem/coefficient.hpp"
#include "fem/errorindicator.hpp"
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/output_functionals.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/rap.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"

namespace palace
{

using namespace Catch::Matchers;

namespace
{

constexpr int crack_attr = 7;
constexpr double x_crack = 0.5;
constexpr double x_sheet = 0.25;  // Second interior boundary, for tests with two of them

// Unit cube mesh of hexahedra for x < x_crack and pyramids for x > x_crack (each hexahedron
// split into six pyramids with their apex at its center).
mfem::Mesh MakeHexPyramidMesh(int n)
{
  auto hex = mfem::Mesh::MakeCartesian3D(n, n, n, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0);
  int num_pyr_hex = 0;
  for (int e = 0; e < hex.GetNE(); e++)
  {
    mfem::Vector c(3);
    hex.GetElementCenter(e, c);
    num_pyr_hex += (c(0) > x_crack);
  }
  mfem::Mesh mesh(3, hex.GetNV() + num_pyr_hex, hex.GetNE() + 5 * num_pyr_hex,
                  hex.GetNBE());
  for (int v = 0; v < hex.GetNV(); v++)
  {
    mesh.AddVertex(hex.GetVertex(v));
  }
  for (int e = 0; e < hex.GetNE(); e++)
  {
    mfem::Vector c(3);
    hex.GetElementCenter(e, c);
    const auto &el = *hex.GetElement(e);
    if (c(0) < x_crack)
    {
      mesh.AddElement(el.Duplicate(&mesh));
      continue;
    }
    const int apex = mesh.AddVertex(c.GetData());
    const int *v = el.GetVertices();
    for (int f = 0; f < 6; f++)
    {
      // The hexahedron faces are oriented outward, the pyramid bases toward the apex.
      const int *fv = mfem::Geometry::Constants<mfem::Geometry::CUBE>::FaceVert[f];
      const int pv[5] = {v[fv[0]], v[fv[3]], v[fv[2]], v[fv[1]], apex};
      mesh.AddElement(new mfem::Pyramid(pv, 1));
    }
  }
  for (int be = 0; be < hex.GetNBE(); be++)
  {
    mesh.AddBdrElement(hex.GetBdrElement(be)->Duplicate(&mesh));
  }
  mesh.FinalizeTopology();
  mesh.Finalize(false, true);
  return mesh;
}

// Unit cube (or square) mesh with interior boundary elements (attribute crack_attr) on the
// plane x = x_crack (or on the given planes), for y <= y_max (the interior boundary has a
// free edge for y_max < 1). For type mfem::Element::PYRAMID, the mesh is that of
// MakeHexPyramidMesh.
mfem::Mesh MakeCrackedCubeMesh(int n, mfem::Element::Type type, double y_max = 1.0,
                               const std::vector<double> &x_planes = {x_crack})
{
  const bool is_2d =
      (type == mfem::Element::TRIANGLE || type == mfem::Element::QUADRILATERAL);
  auto base = is_2d ? mfem::Mesh::MakeCartesian2D(n, n, type, false, 1.0, 1.0)
              : (type == mfem::Element::PYRAMID)
                  ? MakeHexPyramidMesh(n)
                  : mfem::Mesh::MakeCartesian3D(n, n, n, type, 1.0, 1.0, 1.0);
  const int dim = base.Dimension();
  std::vector<int> crack_faces;
  mfem::Array<int> fv;
  for (int f = 0; f < base.GetNumFaces(); f++)
  {
    int e1, e2;
    base.GetFaceElements(f, &e1, &e2);
    if (e2 < 0)
    {
      continue;
    }
    base.GetFaceVertices(f, fv);
    for (auto x_plane : x_planes)
    {
      bool on_plane = true;
      for (auto v : fv)
      {
        on_plane = on_plane && (std::abs(base.GetVertex(v)[0] - x_plane) < 1.0e-12) &&
                   (base.GetVertex(v)[1] < y_max + 1.0e-12);
      }
      if (on_plane)
      {
        crack_faces.push_back(f);
      }
    }
  }
  mfem::Mesh mesh(dim, base.GetNV(), base.GetNE(),
                  base.GetNBE() + static_cast<int>(crack_faces.size()));
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
  for (auto f : crack_faces)
  {
    auto *el = base.GetFace(f)->Duplicate(&mesh);
    el->SetAttribute(crack_attr);
    mesh.AddBdrElement(el);
  }
  mesh.FinalizeTopology();
  mesh.Finalize();
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    mesh.SetAttribute(e, 1);
  }
  mesh.SetAttributes();
  return mesh;
}

// Air bridge standing on a ground plane, as interior boundaries (attribute crack_attr) in
// a unit cube mesh: the ground plane z = 0.25, and a bridge over the strip
// 0.25 < x < 0.75, made of the legs x = 0.25 and x = 0.75 for 0.25 < z < z_top and of the
// top z = z_top, for y_lo < y < y_hi (open at both ends, or a tunnel through the domain for
// y_lo = 0 and y_hi = 1).
mfem::Mesh MakeAirBridgeMesh(int n, mfem::Element::Type type, double y_lo = 0.25,
                             double y_hi = 0.75, double z_top = 0.625)
{
  constexpr double z_ground = 0.25, lo = 0.25, hi = 0.75, eps = 1.0e-12;
  auto base = mfem::Mesh::MakeCartesian3D(n, n, n, type, 1.0, 1.0, 1.0);
  auto In = [](double t, double a, double b) { return t > a - eps && t < b + eps; };
  auto Near = [](double t, double a) { return std::abs(t - a) < eps; };
  const std::array<std::function<bool(const double *)>, 4> pieces = {
      [&](const double *x) { return Near(x[2], z_ground); }, [&](const double *x)
      { return Near(x[2], z_top) && In(x[0], lo, hi) && In(x[1], y_lo, y_hi); },
      [&](const double *x)
      { return Near(x[0], lo) && In(x[2], z_ground, z_top) && In(x[1], y_lo, y_hi); },
      [&](const double *x)
      { return Near(x[0], hi) && In(x[2], z_ground, z_top) && In(x[1], y_lo, y_hi); }};
  std::vector<int> sheet_faces;
  mfem::Array<int> fv;
  for (int f = 0; f < base.GetNumFaces(); f++)
  {
    int e1, e2;
    base.GetFaceElements(f, &e1, &e2);
    if (e2 < 0)
    {
      continue;
    }
    base.GetFaceVertices(f, fv);
    for (const auto &on_piece : pieces)
    {
      if (std::all_of(fv.begin(), fv.end(),
                      [&](int v) { return on_piece(base.GetVertex(v)); }))
      {
        sheet_faces.push_back(f);
        break;
      }
    }
  }
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
    el->SetAttribute(crack_attr);
    mesh.AddBdrElement(el);
  }
  mesh.FinalizeTopology();
  mesh.Finalize();
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    mesh.SetAttribute(e, 1);
  }
  mesh.SetAttributes();
  return mesh;
}

double ElementCentroidX(const mfem::Mesh &mesh, const mfem::Element &el)
{
  double xc = 0.0;
  for (int k = 0; k < el.GetNVertices(); k++)
  {
    xc += mesh.GetVertex(el.GetVertices()[k])[0];
  }
  return xc / el.GetNVertices();
}

// The same mesh, cut along the interior boundary: the vertices on the plane (except for
// those on the free edge y = y_max of the interior boundary) are duplicated for the
// elements (and boundary elements) with x > x_crack. The element ordering is unchanged.
mfem::Mesh CutMesh(const mfem::Mesh &mesh, double y_max = 1.0)
{
  std::vector<int> copy(mesh.GetNV(), -1);
  int nv = mesh.GetNV();
  for (int v = 0; v < mesh.GetNV(); v++)
  {
    if (std::abs(mesh.GetVertex(v)[0] - x_crack) < 1.0e-12 &&
        (y_max >= 1.0 || mesh.GetVertex(v)[1] < y_max - 1.0e-12))
    {
      copy[v] = nv++;
    }
  }
  int num_crack_bdr = 0;
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    num_crack_bdr += (mesh.GetBdrAttribute(be) == crack_attr);
  }
  mfem::Mesh cut(mesh.Dimension(), nv, mesh.GetNE(), mesh.GetNBE() + num_crack_bdr);
  for (int v = 0; v < mesh.GetNV(); v++)
  {
    cut.AddVertex(mesh.GetVertex(v));
  }
  for (int v = 0; v < mesh.GetNV(); v++)
  {
    if (copy[v] >= 0)
    {
      cut.AddVertex(mesh.GetVertex(v));
    }
  }
  auto AddRemapped = [&](const mfem::Element &el, bool bdr)
  {
    auto *new_el = el.Duplicate(&cut);
    int *verts = new_el->GetVertices();
    for (int k = 0; k < new_el->GetNVertices(); k++)
    {
      if (copy[verts[k]] >= 0)
      {
        verts[k] = copy[verts[k]];
      }
    }
    bdr ? cut.AddBdrElement(new_el) : cut.AddElement(new_el);
  };
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    const auto &el = *mesh.GetElement(e);
    if (ElementCentroidX(mesh, el) > x_crack)
    {
      AddRemapped(el, false);
    }
    else
    {
      cut.AddElement(el.Duplicate(&cut));
    }
  }
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    const auto &el = *mesh.GetBdrElement(be);
    if (mesh.GetBdrAttribute(be) == crack_attr)
    {
      cut.AddBdrElement(el.Duplicate(&cut));  // Side x < x_crack
      AddRemapped(el, true);                  // Side x > x_crack
    }
    else if (ElementCentroidX(mesh, el) > x_crack)
    {
      AddRemapped(el, true);
    }
    else
    {
      cut.AddBdrElement(el.Duplicate(&cut));
    }
  }
  cut.FinalizeTopology();
  cut.Finalize();
  cut.SetAttributes();
  return cut;
}

// The same mesh with the elements selected by the predicate on their center numbered
// first.
template <typename Predicate>
mfem::Mesh ReorderElements(const mfem::Mesh &mesh, Predicate &&first)
{
  mfem::Mesh out(mesh.Dimension(), mesh.GetNV(), mesh.GetNE(), mesh.GetNBE());
  for (int v = 0; v < mesh.GetNV(); v++)
  {
    out.AddVertex(mesh.GetVertex(v));
  }
  mfem::Vector c(mesh.SpaceDimension());
  for (bool pass : {true, false})
  {
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      const_cast<mfem::Mesh &>(mesh).GetElementCenter(e, c);
      if (first(c) == pass)
      {
        out.AddElement(mesh.GetElement(e)->Duplicate(&out));
      }
    }
  }
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    out.AddBdrElement(mesh.GetBdrElement(be)->Duplicate(&out));
  }
  out.FinalizeTopology();
  out.Finalize();
  out.SetAttributes();
  return out;
}

// Partition such that process boundaries both coincide with the interior boundary and
// cross it.
std::vector<int> Partition(const mfem::Mesh &mesh, int num_procs)
{
  std::vector<int> part(mesh.GetNE());
  mfem::Vector c(mesh.SpaceDimension());
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    const_cast<mfem::Mesh &>(mesh).GetElementCenter(e, c);
    part[e] =
        ((c(0) > x_crack) + 2 * (c(1) > 0.5) + (c.Size() > 2 && c(2) > 0.5)) % num_procs;
  }
  return part;
}

Mesh MakeParMesh(MPI_Comm comm, mfem::Mesh &smesh, std::vector<int> &part)
{
  smesh.EnsureNodes();
  return Mesh(std::make_unique<mfem::ParMesh>(comm, smesh, part.data()));
}

// True DOFs of the interpolant of a vector field.
Vector Interpolate(FiniteElementSpace &fespace,
                   void (*F)(const mfem::Vector &, mfem::Vector &))
{
  mfem::VectorFunctionCoefficient coeff(fespace.SpaceDimension(), F);
  mfem::ParGridFunction gf(&fespace.Get());
  gf.ProjectCoefficient(coeff);
  Vector x(fespace.GetTrueVSize());
  gf.ParallelProject(x);
  return x;
}

// True DOFs of the interpolant of a scalar field.
Vector Interpolate(FiniteElementSpace &fespace, double (*F)(const mfem::Vector &))
{
  mfem::FunctionCoefficient coeff(F);
  mfem::ParGridFunction gf(&fespace.Get());
  gf.ProjectCoefficient(coeff);
  Vector x(fespace.GetTrueVSize());
  gf.ParallelProject(x);
  return x;
}

double Jump(double x)
{
  return (x > x_crack) ? 1.0 : 0.0;
}

// 2D fields: E with a normal jump, scalar B = ∇ × E with an arbitrary jump.
void SmoothE2D(const mfem::Vector &x, mfem::Vector &E)
{
  E(0) = 1.0 + x(0) * x(1) + 2.0 * Jump(x(0));
  E(1) = std::sin(2.0 * x(1)) + x(0);
}

double SmoothB2D(const mfem::Vector &x)
{
  return 1.0 + x(0) * x(1) * x(1) - 3.0 * Jump(x(0)) * x(1);
}

void JumpE2D(const mfem::Vector &x, mfem::Vector &E)
{
  E(0) = 1.0 + Jump(x(0));
  E(1) = 0.5;
}

double JumpB2D(const mfem::Vector &x)
{
  return 1.0 - 2.0 * Jump(x(0));
}

// Fields with a smooth part and a jump across the interior boundary in the components which
// the respective spaces allow to be discontinuous (normal for ND, tangential for RT).
void SmoothE(const mfem::Vector &x, mfem::Vector &E)
{
  E(0) = 1.0 + x(0) * x(1) + 2.0 * Jump(x(0));
  E(1) = std::sin(2.0 * x(1)) * x(2);
  E(2) = x(0) * x(2) + x(1);
}

void SmoothB(const mfem::Vector &x, mfem::Vector &B)
{
  B(0) = x(1) * x(2) + 1.0;
  B(1) = x(0) + std::cos(x(2)) - 3.0 * Jump(x(0));
  B(2) = x(0) * x(1) + 2.0 * Jump(x(0)) * x(1);
}

// Piecewise constant fields with a jump, recovered exactly by broken spaces.
void JumpE(const mfem::Vector &x, mfem::Vector &E)
{
  E(0) = 1.0 + Jump(x(0));
  E(1) = 0.5;
  E(2) = -0.25;
}

void JumpB(const mfem::Vector &x, mfem::Vector &B)
{
  B(0) = 0.75;
  B(1) = 1.0 - 2.0 * Jump(x(0));
  B(2) = 0.5 * Jump(x(0));
}

// Piecewise constant fields with jumps across the two interior boundaries x = x_sheet and
// x = x_crack (different values in the three regions).
double Jump2(double x)
{
  return ((x > x_sheet) ? 1.0 : 0.0) + 2.0 * Jump(x);
}

void JumpE2(const mfem::Vector &x, mfem::Vector &E)
{
  E(0) = 1.0 + Jump2(x(0));
  E(1) = 0.5;
  E(2) = -0.25;
}

void JumpB2(const mfem::Vector &x, mfem::Vector &B)
{
  B(0) = 0.75;
  B(1) = 1.0 - 2.0 * Jump2(x(0));
  B(2) = 0.5 * Jump2(x(0));
}

// Piecewise constant fields for the tunnel of MakeAirBridgeMesh (y_lo = 0, y_hi = 1, and
// z_top = 0.75): the tangential component B_y (and so H_y) differs in the three regions
// below the ground plane, in the tunnel, and above, with junctions of the three regions at
// the footprints of the legs. The normal component D_z jumps across the ground plane.
int TunnelRegion(const mfem::Vector &x)
{
  if (x(2) < 0.25)
  {
    return 0;
  }
  return (x(0) > 0.25 && x(0) < 0.75 && x(2) < 0.75) ? 1 : 2;
}

void JumpETunnel(const mfem::Vector &x, mfem::Vector &E)
{
  E(0) = 0.5;
  E(1) = -0.25;
  E(2) = (TunnelRegion(x) == 0) ? 1.0 : 2.0;
}

void JumpBTunnel(const mfem::Vector &x, mfem::Vector &B)
{
  constexpr double By[3] = {1.0, -0.5, 2.0};
  B(0) = 0.75;
  B(1) = By[TunnelRegion(x)];
  B(2) = 0.5;
}

// The parallel runs of the estimator tests cover a subset of the orders and element types
// of the serial runs: their parallel aspects (interior boundaries on process boundaries,
// shared and constrained entities, rebalancing and the reconstruction of the sides) do not
// depend on them, and the cost of the parallel runs is dominated by the communication in
// the iterative solves.
template <typename T>
std::vector<T> SerialOrParallel(MPI_Comm comm, std::vector<T> serial,
                                std::vector<T> parallel)
{
  return (Mpi::Size(comm) > 1) ? parallel : serial;
}

struct EstimatorSetup
{
  Mesh mesh;
  MaterialOperator mat_op;
  mfem::ND_FECollection nd_fec;
  mfem::RT_FECollection rt_fec;
  FiniteElementSpaceHierarchy nd_fespaces, rt_fespaces;

  EstimatorSetup(MPI_Comm comm, mfem::Mesh &smesh, std::vector<int> &part, int order,
                 const config::MaterialData &material,
                 const config::PeriodicBoundaryData &periodic)
    : mesh(MakeParMesh(comm, smesh, part)),
      mat_op({material}, periodic, ProblemType::EIGENMODE, mesh), nd_fec(order, 3),
      rt_fec(order - 1, 3),
      nd_fespaces(std::make_unique<FiniteElementSpace>(mesh, &nd_fec)),
      rt_fespaces(std::make_unique<FiniteElementSpace>(mesh, &rt_fec))
  {
  }

  // Element-wise squared estimates for the gradient and curl flux estimators (with Jacobi
  // or algebraic multigrid preconditioning for the flux recovery).
  std::pair<Vector, Vector> Estimate(void (*E)(const mfem::Vector &, mfem::Vector &),
                                     void (*B)(const mfem::Vector &, mfem::Vector &),
                                     const std::vector<int> &crack_attr_list,
                                     bool use_mg = false)
  {
    constexpr double tol = 1.0e-14;
    constexpr int max_it = 10000, print = 0;
    auto &nd_fespace = nd_fespaces.GetFinestFESpace();
    auto &rt_fespace = rt_fespaces.GetFinestFESpace();
    GradFluxErrorEstimator<Vector> grad(mat_op, nd_fespace, rt_fespaces, tol, max_it, print,
                                        use_mg, crack_attr_list);
    CurlFluxErrorEstimator<Vector> curl(mat_op, rt_fespace, nd_fespaces, tol, max_it, print,
                                        use_mg, crack_attr_list);
    ErrorIndicator grad_ind, curl_ind;
    grad.AddErrorIndicator(Interpolate(nd_fespace, E), 0.5, grad_ind);
    curl.AddErrorIndicator(Interpolate(rt_fespace, B), 0.5, curl_ind);
    return {grad_ind.Local(), curl_ind.Local()};
  }
};

struct EstimatorSetup2D
{
  Mesh mesh;
  MaterialOperator mat_op;
  mfem::ND_FECollection nd_fec;
  mfem::RT_FECollection rt_fec;
  mfem::H1_FECollection h1_fec;
  mfem::L2_FECollection l2_fec;
  FiniteElementSpaceHierarchy nd_fespaces, rt_fespaces, h1_fespaces;
  FiniteElementSpace l2_fespace;

  EstimatorSetup2D(MPI_Comm comm, mfem::Mesh &smesh, std::vector<int> &part, int order,
                   const config::MaterialData &material,
                   const config::PeriodicBoundaryData &periodic)
    : mesh(MakeParMesh(comm, smesh, part)),
      mat_op({material}, periodic, ProblemType::EIGENMODE, mesh), nd_fec(order, 2),
      rt_fec(order - 1, 2), h1_fec(order, 2),
      l2_fec(order - 1, 2, mfem::BasisType::GaussLegendre, mfem::FiniteElement::INTEGRAL),
      nd_fespaces(std::make_unique<FiniteElementSpace>(mesh, &nd_fec)),
      rt_fespaces(std::make_unique<FiniteElementSpace>(mesh, &rt_fec)),
      h1_fespaces(std::make_unique<FiniteElementSpace>(mesh, &h1_fec)),
      l2_fespace(mesh, &l2_fec)
  {
  }

  std::pair<Vector, Vector> Estimate(void (*E)(const mfem::Vector &, mfem::Vector &),
                                     double (*B)(const mfem::Vector &),
                                     const std::vector<int> &crack_attr_list)
  {
    constexpr double tol = 1.0e-14;
    constexpr int max_it = 10000, print = 0;
    constexpr bool use_mg = false;
    auto &nd_fespace = nd_fespaces.GetFinestFESpace();
    GradFluxErrorEstimator<Vector> grad(mat_op, nd_fespace, rt_fespaces, tol, max_it, print,
                                        use_mg, crack_attr_list);
    CurlFluxErrorEstimator<Vector> curl(mat_op, l2_fespace, h1_fespaces, tol, max_it, print,
                                        use_mg, crack_attr_list);
    ErrorIndicator grad_ind, curl_ind;
    grad.AddErrorIndicator(Interpolate(nd_fespace, E), 0.5, grad_ind);
    curl.AddErrorIndicator(Interpolate(l2_fespace, B), 0.5, curl_ind);
    return {grad_ind.Local(), curl_ind.Local()};
  }
};

// Compare element-wise estimates, relative to the largest estimate (with an absolute floor
// for estimates which vanish up to round-off).
void CompareEstimates(MPI_Comm comm, const Vector &actual, const Vector &expected,
                      const char *what)
{
  REQUIRE(actual.Size() == expected.Size());
  double max_expected = expected.Size() ? expected.Max() : 0.0;
  Mpi::GlobalMax(1, &max_expected, comm);
  const auto *ha = actual.HostRead();
  const auto *he = expected.HostRead();
  double max_diff = 0.0;
  for (int i = 0; i < actual.Size(); i++)
  {
    max_diff = std::max(max_diff, std::abs(ha[i] - he[i]));
  }
  Mpi::GlobalMax(1, &max_diff, comm);
  INFO(what << ": max. difference " << max_diff << " (max. estimate " << max_expected
            << ")");
  CHECK(max_diff <= 1.0e-8 * max_expected + 1.0e-14);
}

double GlobalMax(MPI_Comm comm, const Vector &v)
{
  double max = v.Size() ? v.Max() : 0.0;
  Mpi::GlobalMax(1, &max, comm);
  return max;
}

bool HasHangingEntities(const mfem::ParMesh &mesh)
{
  int hanging = mesh.Nonconforming() &&
                (mesh.ncmesh->GetEdgeList().slaves.Size() ||
                 (mesh.Dimension() == 3 && mesh.ncmesh->GetFaceList().slaves.Size()));
  Mpi::GlobalMax(1, &hanging, mesh.GetComm());
  return hanging;
}

// Refine the elements selected by the predicate on their center, nonconformingly, as in an
// adaptive loop (without the final mesh update, to allow for rebalancing first).
template <typename Predicate>
void RefineWhere(Mesh &mesh, Predicate &&pred)
{
  auto &pmesh = mesh.Get();
  mfem::Array<int> marked;
  mfem::Vector c(pmesh.SpaceDimension());
  for (int e = 0; e < pmesh.GetNE(); e++)
  {
    pmesh.GetElementCenter(e, c);
    if (pred(c))
    {
      marked.Append(e);
    }
  }
  pmesh.GeneralRefinement(marked, 1, 0);
}

}  // namespace

TEST_CASE("Broken space ignores absent boundary attributes",
          "[brokenspace][Serial][Parallel]")
{
  const auto comm = MPI_COMM_WORLD;
  const auto type = GENERATE(mfem::Element::TRIANGLE, mfem::Element::QUADRILATERAL);
  const bool has_break = GENERATE(false, true);
  auto smesh = MakeCrackedCubeMesh(4, type);
  auto part = Partition(smesh, Mpi::Size(comm));
  auto mesh = MakeParMesh(comm, smesh, part);
  // Test absent attributes both below and above the largest present attribute.
  REQUIRE(mesh.Get().bdr_attributes.Find(6) < 0);
  REQUIRE(mesh.Get().bdr_attributes.Find(9) < 0);
  mfem::H1_FECollection fec(2, 2);
  const std::vector<int> attributes =
      has_break ? std::vector<int>{crack_attr} : std::vector<int>{};
  FiniteElementSpace reference(mesh, &fec, attributes);
  auto extra_attributes = attributes;
  extra_attributes.insert(extra_attributes.end(), {6, 9});
  FiniteElementSpace extra(mesh, &fec, extra_attributes);
  CHECK(extra.GetVSize() == reference.GetVSize());
  CHECK(extra.GetTrueVSize() == reference.GetTrueVSize());
  CHECK(extra.GlobalTrueVSize() == reference.GlobalTrueVSize());

  mfem::Array<int> dofs, reference_dofs;
  mfem::DofTransformation dt;
  int wrong = 0;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    extra.GetElementDofs(e, dofs, dt);
    reference.GetElementDofs(e, reference_dofs, dt);
    wrong +=
        !std::equal(dofs.begin(), dofs.end(), reference_dofs.begin(), reference_dofs.end());
  }
  Mpi::GlobalSum(1, &wrong, comm);
  CHECK(wrong == 0);

  Vector x(reference.GetTrueVSize()), expected(reference.GetVSize()),
      actual(extra.GetVSize());
  for (int t = 0; t < x.Size(); t++)
  {
    x(t) = 1.0 + Mpi::Rank(comm) + 0.01 * t;
  }
  reference.GetProlongationMatrix()->Mult(x, expected);
  extra.GetProlongationMatrix()->Mult(x, actual);
  actual -= expected;
  CHECK(GlobalMax(comm, actual) == 0.0);
}

TEST_CASE("Broken space prolongation", "[brokenspace][Serial][Parallel][GPU]")
{
  const auto comm = MPI_COMM_WORLD;
  const int order = GENERATE(1, 2);
  const bool nonconforming = GENERATE(false, true);
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  fem::DefaultIntegrationOrder::p_trial = order;
  auto smesh = MakeCrackedCubeMesh(4, type);
  if (nonconforming)
  {
    smesh.EnsureNCMesh(true);
  }
  auto part = Partition(smesh, Mpi::Size(comm));
  auto mesh = MakeParMesh(comm, smesh, part);
  const std::vector<int> attr_list = {crack_attr};
  if (nonconforming)
  {
    // Refinement ending on the interior boundary, with hanging entities on it, for which
    // the prolongation is a HypreParMatrix (the sides are computed before the refinement,
    // as by an estimator in an adaptive loop).
    RefineWhere(mesh, [](const mfem::Vector &c)
                { return c(0) > x_crack && c(0) < x_crack + 0.25 && c(1) < 0.5; });
    mesh.Update();
    REQUIRE(HasHangingEntities(mesh.Get()));
  }
  mfem::ND_FECollection nd_fec(order, 3);
  mfem::RT_FECollection rt_fec(order - 1, 3);
  mfem::H1_FECollection h1_fec(order, 3);
  for (const mfem::FiniteElementCollection *fec :
       std::vector<const mfem::FiniteElementCollection *>{&nd_fec, &rt_fec, &h1_fec})
  {
    INFO(fec->Name() << ", hanging " << nonconforming << ", type " << type);
    FiniteElementSpace base_fespace(mesh, fec);
    const auto tsize0 = base_fespace.GlobalTrueVSize();
    FiniteElementSpace fespace(mesh, fec, attr_list);
    REQUIRE(fespace.IsBroken());
    CHECK(fespace.GlobalTrueVSize() > tsize0);

    // An independently owned broken space cannot silently expose its coordinate space.
    CHECK_THROWS_AS(fespace.Get(), mfem::ErrorException);
    CHECK(!base_fespace.IsBroken());
    CHECK(base_fespace.GlobalTrueVSize() == tsize0);
    CHECK(base_fespace.GetProlongationMatrix() != fespace.GetProlongationMatrix());

    // Compare with the global size of the space on the cut mesh.
    if (!nonconforming)
    {
      auto cut = CutMesh(smesh);
      auto cut_mesh = MakeParMesh(comm, cut, part);
      FiniteElementSpace cut_fespace(cut_mesh, fec);
      CHECK(fespace.GlobalTrueVSize() == cut_fespace.GlobalTrueVSize());
    }

    // The vectors are used on the device. The results must not depend on where the output
    // vectors are valid before the products (only on the host, or only on the device with
    // stale contents).
    const Operator &P = *fespace.GetProlongationMatrix();
    const auto *bP = dynamic_cast<const mfem::HypreParMatrix *>(&P);
    REQUIRE(bP);
    CHECK(bP == fespace.Dof_TrueDof_Matrix());
    CHECK(bP->GetGlobalNumRows() == fespace.GlobalVSize());
    CHECK(bP->GetGlobalNumCols() == fespace.GlobalTrueVSize());
    if (nonconforming)
    {
      CHECK(dynamic_cast<const mfem::HypreParMatrix *>(
                base_fespace.GetProlongationMatrix()) != nullptr);
    }
    std::mt19937 gen(Mpi::Rank(comm) + 1);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    Vector x(P.Width()), y(P.Height()), ones(P.Height());
    std::array<Vector, 2> Px, Pty, absPt1;
    for (int k = 0; k < 2; k++)
    {
      Px[k].SetSize(P.Height());
      Pty[k].SetSize(P.Width());
      absPt1[k].SetSize(P.Width());
    }
    for (Vector *v :
         {&x, &y, &ones, &Px[0], &Px[1], &Pty[0], &Pty[1], &absPt1[0], &absPt1[1]})
    {
      v->UseDevice(true);
    }
    for (Vector *v : {&x, &y})
    {
      auto *h_v = v->HostWrite();
      for (int i = 0; i < v->Size(); i++)
      {
        h_v[i] = dist(gen);
      }
    }
    ones = 1.0;
    Px[1] = 1.0e30;
    Pty[1] = 1.0e30;
    absPt1[1] = 1.0e30;
    for (int k = 0; k < 2; k++)
    {
      P.Mult(x, Px[k]);
      P.MultTranspose(y, Pty[k]);
      bP->AbsMultTranspose(1.0, ones, 0.0, absPt1[k]);
    }
    // Reductions on the host, after synchronization (the debug device does not support
    // device reductions).
    auto MaxDiff = [&](const Vector &a, const Vector &b)
    {
      const auto *A = a.HostRead();
      const auto *B = b.HostRead();
      double max = 0.0;
      for (int i = 0; i < a.Size(); i++)
      {
        max = std::max(max, std::abs(A[i] - B[i]));
      }
      Mpi::GlobalMax(1, &max, comm);
      return max;
    };
    auto Dot = [](const Vector &a, const Vector &b)
    {
      const auto *A = a.HostRead();
      const auto *B = b.HostRead();
      double dot = 0.0;
      for (int i = 0; i < a.Size(); i++)
      {
        dot += A[i] * B[i];
      }
      return dot;
    };
    CHECK(MaxDiff(Px[0], Px[1]) <= 1.0e-12);
    CHECK(MaxDiff(Pty[0], Pty[1]) <= 1.0e-12);
    CHECK(MaxDiff(absPt1[0], absPt1[1]) == 0.0);

    // Adjointness of the prolongation and its transpose. The L-vector inner product is not
    // a global inner product (shared L-DOFs are counted on each process), but
    // y · P x = Pᵀ y · x holds summed over processes.
    double dots[2] = {Dot(y, Px[0]), Dot(Pty[0], x)};
    Mpi::GlobalSum(2, dots, comm);
    CHECK_THAT(dots[0], WithinRel(dots[1], 1.0e-12));

    // Accumulated column weight is positive for Jacobi; individual rows need not be unit.
    double min = 1.0, max = 0.0;
    {
      const auto *D = absPt1[0].HostRead();
      for (int i = 0; i < absPt1[0].Size(); i++)
      {
        min = std::min(min, D[i]);
        max = std::max(max, D[i]);
      }
    }
    Mpi::GlobalMin(1, &min, comm);
    Mpi::GlobalMax(1, &max, comm);
    CHECK(min >= 1.0 - 1.0e-12);

    // y = a |P|ᵀ x + b y.
    bP->AbsMultTranspose(2.0, ones, 0.5, absPt1[1]);
    absPt1[0] *= 2.5;
    CHECK(MaxDiff(absPt1[0], absPt1[1]) <= 2.5e-12 * max);

    // Both restriction and DG reverse-copy must recover arbitrary broken fields, not
    // just fields equal on both sides. Includes copies owned only through remote rows.
    Vector recovered(P.Width());
    const auto &R = *fespace.GetRestrictionMatrix();
    R.Mult(Px[0], recovered);
    CHECK(MaxDiff(x, recovered) < 1.0e-12);
    // R A P must use the same local sparse restriction in apply and full assembly.
    auto identity = std::make_unique<mfem::SparseMatrix>(P.Height());
    for (int l = 0; l < P.Height(); l++)
    {
      identity->Set(l, l, 1.0);
    }
    identity->Finalize();
    ParOperator identity_true(std::move(identity), fespace, fespace, true);
    identity_true.Mult(x, recovered);
    CHECK(MaxDiff(x, recovered) < 1.0e-12);
    identity_true.ParallelAssemble().Mult(x, recovered);
    CHECK(MaxDiff(x, recovered) < 1.0e-12);
    Vector Ry(P.Width()), Rtx(P.Height());
    R.Mult(y, Ry);
    R.MultTranspose(x, Rtx);
    double rdots[2] = {Dot(y, Rtx), Dot(Ry, x)};
    Mpi::GlobalSum(2, rdots, comm);
    CHECK_THAT(rdots[0], WithinRel(rdots[1], 1.0e-12) || WithinAbs(rdots[1], 1.0e-12));
    CHECK_THROWS_AS(GridFunction(fespace), mfem::ErrorException);
    GridFunction dg(fespace, GridFunction::Storage::DG);
    dg.SetFromTrueDofs(x);
    dg.GetTrueDofs(recovered);
    CHECK(MaxDiff(x, recovered) < 1.0e-12);
    ComplexVector cx(x.Size()), cy(x.Size());
    cx.Real() = x;
    cx.Imag() = x;
    cx.Imag() *= -0.7;
    GridFunction complex_dg(fespace, GridFunction::Storage::DG, true);
    complex_dg.SetFromTrueDofs(cx);
    complex_dg.GetTrueDofs(cy);
    CHECK(MaxDiff(cx.Real(), cy.Real()) < 1.0e-12);
    CHECK(MaxDiff(cx.Imag(), cy.Imag()) < 1.0e-12);

    // A round trip alone could cancel a missing DOF transformation in both directions.
    // Check physical values against the source element's native shapes as well.
    mfem::Array<int> dofs;
    mfem::DofTransformation dt;
    Vector native, actual, expected;
    mfem::DenseMatrix shape;
    double error = 0.0;
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      const auto *fe = fespace.GetFE(e);
      const auto &ip = mfem::Geometries.GetCenter(fe->GetGeomType());
      auto *T = mesh.Get().GetElementTransformation(e);
      T->SetIntPoint(&ip);
      fespace.GetElementDofs(e, dofs, dt);
      Px[0].GetSubVector(dofs, native);
      if (!dt.IsIdentity())
      {
        dt.InvTransformPrimal(native);
      }
      if (fe->GetRangeType() == mfem::FiniteElement::VECTOR)
      {
        shape.SetSize(fe->GetDof(), mesh.SpaceDimension());
        fe->CalcVShape(*T, shape);
        expected.SetSize(mesh.SpaceDimension());
        shape.MultTranspose(native, expected);
        dg.Real().GetVectorValue(e, ip, actual);
        actual -= expected;
        error = std::max(error, actual.Normlinf());
      }
      else
      {
        Vector scalar_shape(fe->GetDof());
        fe->CalcPhysShape(*T, scalar_shape);
        error =
            std::max(error, std::abs(dg.Real().GetValue(e, ip) - (scalar_shape * native)));
      }
    }
    Mpi::GlobalMax(1, &error, comm);
    CHECK(error < 1.0e-12);

    // A mutable MFEM consumer must not become a DG-to-source import path.
    dg.Real() = 42.0;
    dg.GetTrueDofs(recovered);
    CHECK(MaxDiff(x, recovered) == 0.0);
    recovered.SetSize(fespace.GetTrueVSize());
    fespace.GetRestrictionMatrix()->Mult(dg.GetLVector(), recovered);
    CHECK(MaxDiff(x, recovered) < 1.0e-12);

    // Exercise native tensor interpolation restrictions and the ND tetrahedral range
    // transformation directly; assembling interpolators through the coordinate space
    // would not exercise these public restrictions.
    for (auto ceed : ceed::internal::GetCeedObjects())
    {
      for (const auto &[geom, data] : mesh.GetCeedGeomFactorData(ceed))
      {
        if (mfem::Geometry::Dimension[geom] != mesh.Dimension())
        {
          continue;
        }
        const int n = fespace.GetFEColl().FiniteElementForGeometry(geom)->GetDof();
        for (bool range : {false, true})
        {
          auto restriction =
              range ? fespace.GetInterpRangeCeedElemRestriction(ceed, geom, data)
                    : fespace.GetInterpCeedElemRestriction(ceed, geom, data);
          CeedVector lv, ev;
          PalaceCeedCall(ceed, CeedElemRestrictionCreateVector(restriction, &lv, &ev));
          // This is a new vector, not an existing borrowed array to repoint.
          PalaceCeedCall(ceed,
                         CeedVectorSetArray(lv, CEED_MEM_HOST, CEED_COPY_VALUES,
                                            const_cast<CeedScalar *>(Px[0].HostRead())));
          PalaceCeedCall(ceed, CeedElemRestrictionApply(restriction, CEED_NOTRANSPOSE, lv,
                                                        ev, CEED_REQUEST_IMMEDIATE));
          const CeedScalar *values;
          PalaceCeedCall(ceed, CeedVectorGetArrayRead(ev, CEED_MEM_HOST, &values));
          double restriction_error = 0.0;
          for (std::size_t i = 0; i < data.indices.size(); i++)
          {
            fespace.GetElementDofs(data.indices[i], dofs, dt);
            Px[0].GetSubVector(dofs, native);
            if (!dt.IsIdentity())
            {
              if (range)
              {
                dt.InvTransformDual(native);
              }
              else
              {
                dt.InvTransformPrimal(native);
              }
            }
            for (int j = 0; j < n; j++)
            {
              restriction_error =
                  std::max(restriction_error, std::abs(values[i * n + j] - native(j)));
            }
          }
          PalaceCeedCall(ceed, CeedVectorRestoreArrayRead(ev, &values));
          PalaceCeedCall(ceed, CeedVectorDestroy(&lv));
          PalaceCeedCall(ceed, CeedVectorDestroy(&ev));
          CHECK(restriction_error < 1.0e-12);
        }
      }
    }
  }
}

TEST_CASE("Broken space interior traces", "[brokenspace][Serial][Parallel]")
{
  const auto comm = MPI_COMM_WORLD;
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  const int order = GENERATE(1, 2);
  const bool hanging = GENERATE(false, true);
  auto smesh = MakeCrackedCubeMesh(2, type);
  smesh.EnsureNCMesh(true);
  auto part = Partition(smesh, Mpi::Size(comm));
  auto mesh = MakeParMesh(comm, smesh, part);
  if (hanging)
  {
    RefineWhere(mesh, [](const mfem::Vector &c) { return c(0) > x_crack; });
    mesh.Update();
  }
  mfem::ND_FECollection fec(order, 3);
  FiniteElementSpace space(mesh, &fec, {crack_attr});
  Vector local(space.GetVSize()), true_values(space.GetTrueVSize()), native, center,
      value(3);
  local = 0.0;
  mfem::Array<int> dofs;
  mfem::DofTransformation dt;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    mesh.Get().GetElementCenter(e, center);
    const double side = center(0) < x_crack ? 1.0 : -2.0;
    value = 0.0;
    value(1) = side;
    value(2) = 0.2 * side;
    mfem::VectorConstantCoefficient coeff(value);
    native.SetSize(space.GetFE(e)->GetDof());
    space.GetFE(e)->Project(coeff, *mesh.Get().GetElementTransformation(e), native);
    space.GetElementDofs(e, dofs, dt);
    if (!dt.IsIdentity())
    {
      dt.TransformPrimal(native);
    }
    local.SetSubVector(dofs, native);
  }
  space.GetRestrictionMatrix()->Mult(local, true_values);
  space.GetProlongationMatrix()->Mult(true_values, local);
  const auto faces = space.GetInteriorBoundaryFaces({9999, crack_attr, crack_attr});
  double error = 0.0, area = 0.0;
  for (auto f : faces)
  {
    const auto info = mesh.Get().GetFaceInformation(f);
    auto *T = info.IsShared() ? mesh.Get().GetSharedFaceTransformationsByLocalIndex(f)
                              : mesh.Get().GetFaceElementTransformations(f);
    mfem::Array<int> first, second;
    mfem::DofTransformation dt1, dt2;
    space.GetInteriorFaceDofs(f, first, dt1, second, dt2);
    const auto *fe1 = space.GetFE(info.element[0].index);
    const auto *fe2 = info.IsShared() ? space.GetFaceNbrFE(info.element[1].index)
                                      : space.GetFE(info.element[1].index);
    const auto &ir = mfem::IntRules.Get(T->GetGeometryType(), 2 * order + 3);
    for (int q = 0; q < ir.GetNPoints(); q++)
    {
      T->SetAllIntPoints(&ir.IntPoint(q));
      area += ir.IntPoint(q).weight * T->Weight();
      auto Check = [&](const mfem::FiniteElement &fe, mfem::ElementTransformation &element,
                       const mfem::Array<int> &read, mfem::DofTransformation &transform)
      {
        local.GetSubVector(read, native);
        if (!transform.IsIdentity())
        {
          transform.InvTransformPrimal(native);
        }
        mfem::DenseMatrix shape(fe.GetDof(), 3);
        fe.CalcVShape(element, shape);
        shape.MultTranspose(native, value);
        element.Transform(mfem::Geometries.GetCenter(fe.GetGeomType()), center);
        const double expected = center(0) < x_crack ? 1.0 : -2.0;
        error = std::max({error, std::abs(value(0)), std::abs(value(1) - expected),
                          std::abs(value(2) - 0.2 * expected)});
      };
      Check(*fe1, *T->Elem1, first, dt1);
      Check(*fe2, *T->Elem2, second, dt2);
    }
  }
  Mpi::GlobalMax(1, &error, comm);
  Mpi::GlobalSum(1, &area, comm);
  CHECK(error < 1.0e-11);
  CHECK_THAT(area, WithinAbs(1.0, 1.0e-12));
}

TEST_CASE("Broken restriction on an edge-only reading rank", "[brokenspace][Parallel]")
{
  const auto comm = MPI_COMM_WORLD;
  if (Mpi::Size(comm) < 3)
  {
    SKIP("Needs three ranks");
  }
  auto smesh = MakeCrackedCubeMesh(2, mfem::Element::TETRAHEDRON);
  smesh.EnsureNCMesh(true);
  mfem::Array<int> refined, vertices;
  Vector center;
  for (int e = 0; e < smesh.GetNE(); e++)
  {
    smesh.GetElementCenter(e, center);
    if (center(0) > x_crack)
    {
      refined.Append(e);
    }
  }
  smesh.GeneralRefinement(refined, 1, 0);
  std::vector<int> part(smesh.GetNE());
  int selected = -1;
  for (int e = 0; e < smesh.GetNE(); e++)
  {
    smesh.GetElementCenter(e, center);
    part[e] = center(0) < x_crack ? 2 : 1;
    if (center(0) < x_crack)
    {
      continue;
    }
    smesh.GetElementVertices(e, vertices);
    std::vector<int> plane;
    for (auto v : vertices)
    {
      if (std::abs(smesh.GetVertex(v)[0] - x_crack) < 1.0e-14)
      {
        plane.push_back(v);
      }
    }
    if (selected < 0 && plane.size() == 2 &&
        std::abs(smesh.GetVertex(plane[0])[1] - 0.25) < 1.0e-14 &&
        std::abs(smesh.GetVertex(plane[1])[1] - 0.25) < 1.0e-14)
    {
      selected = e;
      part[e] = 0;
    }
  }
  REQUIRE(selected >= 0);
  auto mesh = MakeParMesh(comm, smesh, part);
  mfem::H1_FECollection h1_fec(2, 3);
  mfem::ND_FECollection nd_fec(2, 3);
  mfem::RT_FECollection rt_fec(1, 3);
  FiniteElementSpace h1(mesh, &h1_fec, {crack_attr}), nd(mesh, &nd_fec, {crack_attr}),
      rt(mesh, &rt_fec, {crack_attr});
  Vector x(h1.GetTrueVSize()), gradient(nd.GetTrueVSize()), curl(rt.GetTrueVSize());
  for (int i = 0; i < x.Size(); i++)
  {
    x(i) = std::sin(i + 1.7 * Mpi::Rank(comm));
  }
  nd.GetDiscreteInterpolator(h1).Mult(x, gradient);
  rt.GetDiscreteInterpolator(nd).Mult(gradient, curl);
  double error = curl.Normlinf(), scale = gradient.Normlinf();
  Mpi::GlobalMax(1, &error, comm);
  Mpi::GlobalMax(1, &scale, comm);
  CHECK(error < 1.0e-12 * scale);
  Vector input(nd.GetTrueVSize()), local(nd.GetVSize()), recovered(nd.GetTrueVSize());
  for (int i = 0; i < input.Size(); i++)
  {
    input(i) = std::cos(0.8 * i + Mpi::Rank(comm));
  }
  rt.GetDiscreteInterpolator(nd).Mult(input, curl);
  nd.GetProlongationMatrix()->Mult(input, local);
  nd.GetRestrictionMatrix()->Mult(local, recovered);
  recovered -= input;
  error = recovered.Normlinf();
  Mpi::GlobalMax(1, &error, comm);
  CHECK(error < 1.0e-12 * scale);

  Mesh fine(std::make_unique<mfem::ParMesh>(mesh.Get()));
  fine.Get().UniformRefinement();
  FiniteElementSpaceHierarchy nh, rh;
  nh.AddLevel(
      std::make_unique<FiniteElementSpace>(mesh, &nd_fec, std::vector<int>{crack_attr}));
  nh.AddLevel(
      std::make_unique<FiniteElementSpace>(fine, &nd_fec, std::vector<int>{crack_attr}));
  rh.AddLevel(
      std::make_unique<FiniteElementSpace>(mesh, &rt_fec, std::vector<int>{crack_attr}));
  rh.AddLevel(
      std::make_unique<FiniteElementSpace>(fine, &rt_fec, std::vector<int>{crack_attr}));
  Vector refined_nd(nh.GetFinestFESpace().GetTrueVSize()),
      refined_rt(rh.GetFinestFESpace().GetTrueVSize()), expected(refined_rt.Size());
  nh.GetProlongationAtLevel(0).Mult(input, refined_nd);
  rh.GetFinestFESpace()
      .GetDiscreteInterpolator(nh.GetFinestFESpace())
      .Mult(refined_nd, refined_rt);
  rh.GetProlongationAtLevel(0).Mult(curl, expected);
  refined_rt -= expected;
  error = refined_rt.Normlinf();
  Mpi::GlobalMax(1, &error, comm);
  CHECK(error < 1.0e-11 * scale);
}

TEST_CASE("Broken space hanging sheet version counts", "[brokenspace][Serial][Parallel]")
{
  const auto comm = MPI_COMM_WORLD;
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  const int order = GENERATE(1, 2);
  auto smesh = MakeCrackedCubeMesh(2, type);
  smesh.EnsureNCMesh(true);
  auto part = Partition(smesh, Mpi::Size(comm));
  auto mesh = MakeParMesh(comm, smesh, part);
  mfem::H1_FECollection h1(order, 3);
  mfem::ND_FECollection nd(order, 3);
  mfem::RT_FECollection rt(order - 1, 3);
  std::vector<const mfem::FiniteElementCollection *> collections = {&h1, &nd, &rt};
  std::vector<HYPRE_BigInt> copies;
  for (auto *fec : collections)
  {
    FiniteElementSpace ordinary(mesh, fec), broken(mesh, fec, {crack_attr});
    copies.push_back(broken.GlobalTrueVSize() - ordinary.GlobalTrueVSize());
    REQUIRE(copies.back() > 0);
  }
  for (int level = 0; level < 2; level++)
  {
    RefineWhere(mesh, [](const mfem::Vector &c) { return c(0) > x_crack; });
    mesh.Get().Rebalance();
    mesh.Update();
    REQUIRE(HasHangingEntities(mesh.Get()));
    for (std::size_t i = 0; i < collections.size(); i++)
    {
      INFO(collections[i]->Name() << ", level " << level << ", type " << type);
      // Both traces retain the master's trace space. In particular its copies on
      // the fine side are NOT independent copies for each child face.
      FiniteElementSpace ordinary(mesh, collections[i]);
      FiniteElementSpace broken(mesh, collections[i], {crack_attr});
      CHECK(broken.GlobalTrueVSize() - ordinary.GlobalTrueVSize() == copies[i]);
    }
  }
}

TEST_CASE("Broken space boundary DOFs", "[brokenspace][Serial][Parallel]")
{
  // A boundary standing on one side of the interior boundary (the sheet y = 0.5 for
  // x > x_crack, attribute 9, which meets the interior boundary along an edge) reads and
  // makes essential the versions of that side only: each of its DOFs in the underlying
  // space is essential in exactly one version of the broken space.
  const auto comm = MPI_COMM_WORLD;
  const int order = GENERATE(1, 2);
  constexpr int sheet_attr = 9;
  auto base = MakeCrackedCubeMesh(4, mfem::Element::TETRAHEDRON);
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
      on_sheet = on_sheet && std::abs(base.GetVertex(v)[1] - 0.5) < 1.0e-12 &&
                 base.GetVertex(v)[0] > x_crack - 1.0e-12;
    }
    if (on_sheet)
    {
      sheet_faces.push_back(f);
    }
  }
  REQUIRE(!sheet_faces.empty());
  mfem::Mesh smesh(3, base.GetNV(), base.GetNE(),
                   base.GetNBE() + static_cast<int>(sheet_faces.size()));
  for (int v = 0; v < base.GetNV(); v++)
  {
    smesh.AddVertex(base.GetVertex(v));
  }
  for (int e = 0; e < base.GetNE(); e++)
  {
    smesh.AddElement(base.GetElement(e)->Duplicate(&smesh));
  }
  for (int be = 0; be < base.GetNBE(); be++)
  {
    smesh.AddBdrElement(base.GetBdrElement(be)->Duplicate(&smesh));
  }
  for (auto f : sheet_faces)
  {
    auto *el = base.GetFace(f)->Duplicate(&smesh);
    el->SetAttribute(sheet_attr);
    smesh.AddBdrElement(el);
  }
  smesh.FinalizeTopology();
  smesh.Finalize();
  smesh.SetAttributes();
  smesh.EnsureNCMesh(true);
  auto part = Partition(smesh, Mpi::Size(comm));
  auto mesh = MakeParMesh(comm, smesh, part);
  mfem::ND_FECollection nd_fec(order, 3);
  FiniteElementSpace nd_fespace(mesh, &nd_fec);
  FiniteElementSpace nd_broken(mesh, &nd_fec, {crack_attr});
  REQUIRE(nd_broken.IsBroken());
  CHECK(nd_broken.GetVDim() == 1);
  mfem::Array<int> marker(mesh.Get().bdr_attributes.Max());
  marker = 0;
  marker[sheet_attr - 1] = 1;

  // The boundary element DOFs are those read by the adjacent element.
  mfem::Array<int> bdofs, edofs;
  mfem::DofTransformation dof_trans;
  int num_bdr = 0;
  for (int be = 0; be < mesh.Get().GetNBE(); be++)
  {
    if (mesh.Get().GetBdrAttribute(be) != sheet_attr)
    {
      continue;
    }
    int e, info;
    mesh.Get().GetBdrElementAdjacentElement(be, e, info);
    nd_broken.GetBdrElementDofs(be, bdofs, dof_trans);
    nd_broken.GetElementDofs(e, edofs, dof_trans);
    for (auto d : bdofs)
    {
      const int l = (d >= 0) ? d : -1 - d;
      bool found = false;
      for (auto ed : edofs)
      {
        found = found || ((ed >= 0) ? ed : -1 - ed) == l;
      }
      CHECK(found);
    }
    num_bdr++;
  }
  Mpi::GlobalSum(1, &num_bdr, comm);
  CHECK(num_bdr == static_cast<int>(sheet_faces.size()));

  // One version per essential DOF of the underlying space, and only the versions of the
  // sheet side: the elements on the other side read no essential L-DOF (the sheet's DOFs on
  // the interior boundary are split, since the interior boundary spans the domain there).
  auto CheckEssential = [&]()
  {
    mfem::Array<int> ess_u, ess_b;
    nd_fespace.GetEssentialTrueDofs(marker, ess_u);
    nd_broken.GetEssentialTrueDofs(marker, ess_b);
    int sizes[2] = {ess_u.Size(), ess_b.Size()};
    Mpi::GlobalSum(2, sizes, comm);
    CHECK(sizes[0] == sizes[1]);
    CHECK(nd_broken.GlobalTrueVSize() > nd_fespace.GlobalTrueVSize());
    Vector ess_t(nd_broken.GetTrueVSize()), ess_l(nd_broken.GetVSize());
    ess_t.UseDevice(false);
    ess_l.UseDevice(false);
    ess_t = 0.0;
    for (auto t : ess_b)
    {
      ess_t(t) = 1.0;
    }
    nd_broken.GetProlongationMatrix()->Mult(ess_t, ess_l);
    int num_wrong = 0;
    for (int e = 0; e < mesh.Get().GetNE(); e++)
    {
      mfem::Vector c;
      mesh.Get().GetElementCenter(e, c);
      if (c(0) > x_crack)
      {
        continue;
      }
      nd_broken.GetElementDofs(e, edofs, dof_trans);
      for (auto d : edofs)
      {
        num_wrong += (ess_l((d >= 0) ? d : -1 - d) != 0.0);
      }
    }
    Mpi::GlobalSum(1, &num_wrong, comm);
    CHECK(num_wrong == 0);
  };
  CheckEssential();

  // The actual libCEED boundary operator, not only the DOF accessor, must read that
  // side. Compare with independently assembled MFEM element matrices on the L-vector.
  fem::DefaultIntegrationOrder::p_trial = order;
  auto CheckBoundary = [&]()
  {
    BilinearForm boundary(nd_broken);
    boundary.AddBoundaryIntegrator<VectorFEMassIntegrator>();
    auto op = boundary.PartialAssemble();
    Vector local(nd_broken.GetVSize()), result(local.Size()), reference(local.Size());
    for (int i = 0; i < local.Size(); i++)
    {
      local(i) = std::sin(1.7 * (i + 1 + Mpi::Rank(comm)));
    }
    reference = 0.0;
    mfem::VectorFEMassIntegrator mass;
    mfem::DenseMatrix matrix;
    Vector element;
    for (int be = 0; be < mesh.Get().GetNBE(); be++)
    {
      auto *T = mesh.Get().GetBdrElementTransformation(be);
      const auto &ir =
          mfem::IntRules.Get(T->GetGeometryType(), fem::DefaultIntegrationOrder::Get(*T));
      mass.SetIntRule(&ir);
      mass.AssembleElementMatrix(*nd_broken.GetBE(be), *T, matrix);
      nd_broken.GetBdrElementDofs(be, bdofs, dof_trans);
      local.GetSubVector(bdofs, element);
      if (!dof_trans.IsIdentity())
      {
        dof_trans.InvTransformPrimal(element);
      }
      Vector value(element.Size());
      matrix.Mult(element, value);
      if (!dof_trans.IsIdentity())
      {
        dof_trans.TransformDual(value);
      }
      reference.AddElementVector(bdofs, value);
    }
    op->Mult(local, result);
    result -= reference;
    double error = result.Normlinf();
    Mpi::GlobalMax(1, &error, comm);
    CHECK(error < 1.0e-12);
  };
  CheckBoundary();

  // Both independent break lists are recomputed on the current mesh after each update.
  // Adding the standing sheet gives a three-sided junction, allowed for a general space.
  FiniteElementSpace more_breaks(mesh, &nd_fec, std::vector<int>{sheet_attr, crack_attr});
  CHECK(more_breaks.GlobalTrueVSize() > nd_broken.GlobalTrueVSize());
  auto CheckTwoLists = [&]()
  {
    FiniteElementSpace rebuilt(mesh, &nd_fec, std::vector<int>{crack_attr});
    FiniteElementSpace rebuilt_more(mesh, &nd_fec,
                                    std::vector<int>{crack_attr, sheet_attr, sheet_attr});
    CHECK(rebuilt.GlobalTrueVSize() == nd_broken.GlobalTrueVSize());
    CHECK(rebuilt_more.GlobalTrueVSize() == more_breaks.GlobalTrueVSize());
    CHECK(more_breaks.GetBreakAttributes() == rebuilt_more.GetBreakAttributes());
  };
  CheckTwoLists();
  mfem::H1_FECollection h1_fec(order, 3);
  FiniteElementSpace wrong_aux(mesh, &h1_fec, std::vector<int>{sheet_attr, crack_attr});
  CHECK_THROWS_AS(nd_broken.GetDiscreteInterpolator(wrong_aux), mfem::ErrorException);
  RefineWhere(mesh, [](const mfem::Vector &c) { return c(0) > x_crack && c(1) < 0.5; });
  mesh.Update();
  nd_fespace.Update();
  nd_broken.Update();
  more_breaks.Update();
  REQUIRE(HasHangingEntities(mesh.Get()));
  CheckEssential();
  CheckTwoLists();
  CheckBoundary();
  mesh.Get().Rebalance();
  mesh.Update();
  nd_fespace.Update();
  nd_broken.Update();
  more_breaks.Update();
  CheckEssential();
  CheckTwoLists();
  CheckBoundary();
}

TEST_CASE("Broken space interpolators", "[brokenspace][Serial][Parallel]")
{
  // The discrete gradient, curl, and p-prolongation between spaces broken with the same
  // sides, which are assembled from the interpolators between the underlying spaces (see
  // FiniteElementSpace::BuildBrokenInterpolator), including on a mesh with hanging
  // entities on the interior boundary, where some copies are read only through the
  // constraints.
  const auto comm = MPI_COMM_WORLD;
  const int order = GENERATE(1, 2);
  const bool nonconforming = GENERATE(false, true);
  const bool free_edge = GENERATE(false, true);
  CAPTURE(order, nonconforming, free_edge);
  fem::DefaultIntegrationOrder::p_trial = order;
  auto smesh = MakeCrackedCubeMesh(4, mfem::Element::TETRAHEDRON, free_edge ? 0.5 : 1.0);
  if (nonconforming)
  {
    smesh.EnsureNCMesh(true);
  }
  auto part = Partition(smesh, Mpi::Size(comm));
  auto mesh = MakeParMesh(comm, smesh, part);
  const std::vector<int> attr_list = {crack_attr};
  if (nonconforming)
  {
    RefineWhere(mesh, [](const mfem::Vector &c)
                { return c(0) > x_crack && c(0) < x_crack + 0.25 && c(1) < 0.5; });
    mesh.Update();
    REQUIRE(HasHangingEntities(mesh.Get()));
  }
  mfem::H1_FECollection h1_fec(order, 3);
  mfem::ND_FECollection nd_fec(order, 3), nd_fec_coarse(std::max(order - 1, 1), 3);
  mfem::RT_FECollection rt_fec(order - 1, 3), rt_fec_coarse(std::max(order - 2, 0), 3);
  FiniteElementSpace h1_fespace(mesh, &h1_fec), nd_fespace(mesh, &nd_fec),
      rt_fespace(mesh, &rt_fec);
  FiniteElementSpace h1_broken(mesh, &h1_fec, attr_list),
      nd_broken(mesh, &nd_fec, attr_list), rt_broken(mesh, &rt_fec, attr_list);
  FiniteElementSpaceHierarchy nd_hierarchy, rt_hierarchy;
  nd_hierarchy.AddLevel(
      std::make_unique<FiniteElementSpace>(mesh, &nd_fec_coarse, attr_list));
  nd_hierarchy.AddLevel(std::make_unique<FiniteElementSpace>(mesh, &nd_fec, attr_list));
  rt_hierarchy.AddLevel(
      std::make_unique<FiniteElementSpace>(mesh, &rt_fec_coarse, attr_list));
  rt_hierarchy.AddLevel(std::make_unique<FiniteElementSpace>(mesh, &rt_fec, attr_list));
  const Operator &G = nd_broken.GetDiscreteInterpolator(h1_broken);
  const Operator &C = rt_broken.GetDiscreteInterpolator(nd_broken);
  const Operator &P_nd = nd_hierarchy.GetProlongationAtLevel(0);
  const Operator &P_rt = rt_hierarchy.GetProlongationAtLevel(0);
  const Operator &C_coarse = rt_hierarchy.GetFESpaceAtLevel(0).GetDiscreteInterpolator(
      nd_hierarchy.GetFESpaceAtLevel(0));
  CHECK(dynamic_cast<const mfem::HypreParMatrix *>(&G));
  CHECK(G.Height() == nd_broken.GetTrueVSize());
  CHECK(G.Width() == h1_broken.GetTrueVSize());

  std::mt19937 gen(Mpi::Rank(comm) + 1);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  auto Random = [&](int n)
  {
    Vector v(n);
    for (int i = 0; i < n; i++)
    {
      v(i) = dist(gen);
    }
    return v;
  };
  auto MaxAbs = [&](const Vector &v)
  {
    double max = v.Normlinf();
    Mpi::GlobalMax(1, &max, comm);
    return max;
  };

  // A field of the underlying space, copied to all versions, is interpolated as by the
  // interpolator of the underlying spaces.
  FiniteElementSpaceHierarchy nd_hierarchy_u;
  nd_hierarchy_u.AddLevel(std::make_unique<FiniteElementSpace>(mesh, &nd_fec_coarse));
  nd_hierarchy_u.AddLevel(std::make_unique<FiniteElementSpace>(mesh, &nd_fec));
  auto CheckCopies = [&](const FiniteElementSpace &trial, const FiniteElementSpace &test,
                         const Operator &G_u, const Operator &G_b)
  {
    Vector x = Random(G_u.Width()), y(G_u.Height()), xb(G_b.Width()), yb(G_b.Height()),
           zb(G_b.Height());
    G_u.Mult(x, y);
    trial.CopyFromConforming(x, xb);
    G_b.Mult(xb, yb);
    test.CopyFromConforming(y, zb);
    zb -= yb;
    CHECK(MaxAbs(zb) <= 1.0e-12 * MaxAbs(yb));
  };
  CheckCopies(h1_broken, nd_broken, nd_fespace.GetDiscreteInterpolator(h1_fespace), G);
  CheckCopies(nd_broken, rt_broken, rt_fespace.GetDiscreteInterpolator(nd_fespace), C);
  CheckCopies(nd_hierarchy.GetFESpaceAtLevel(0), nd_broken,
              nd_hierarchy_u.GetProlongationAtLevel(0), P_nd);

  // The de Rham complex of the broken spaces: the curl of a gradient of any broken field
  // (with different values on the two sides) is zero, and the curl of the p-prolongation
  // of a broken field is the p-prolongation of its curl.
  {
    Vector x = Random(G.Width()), g(G.Height()), c(C.Height());
    G.Mult(x, g);
    C.Mult(g, c);
    CHECK(MaxAbs(c) <= 1.0e-12 * MaxAbs(g));
  }
  {
    Vector x = Random(P_nd.Width()), px(P_nd.Height()), cpx(C.Height()),
           cx(C_coarse.Height()), pcx(P_rt.Height());
    P_nd.Mult(x, px);
    C.Mult(px, cpx);
    C_coarse.Mult(x, cx);
    P_rt.Mult(cx, pcx);
    pcx -= cpx;
    CHECK(MaxAbs(pcx) <= 1.0e-12 * MaxAbs(cpx));
  }

  // Two mesh levels: the fine mesh is a uniform refinement of the mesh, with the sides
  // inherited from it, and the h-prolongations of broken hierarchies are assembled from
  // MFEM's transfer operators. The curl of the h-prolongation of a broken field is the
  // h-prolongation of its curl.
  if (!nonconforming)
  {
    Mesh fine_mesh(std::make_unique<mfem::ParMesh>(mesh.Get()));
    fine_mesh.Get().UniformRefinement();
    FiniteElementSpace nd_fine(fine_mesh, &nd_fec), rt_fine(fine_mesh, &rt_fec);
    FiniteElementSpaceHierarchy nd_h, rt_h, nd_h_u;
    nd_h.AddLevel(std::make_unique<FiniteElementSpace>(mesh, &nd_fec, attr_list));
    nd_h.AddLevel(std::make_unique<FiniteElementSpace>(fine_mesh, &nd_fec, attr_list));
    rt_h.AddLevel(std::make_unique<FiniteElementSpace>(mesh, &rt_fec, attr_list));
    rt_h.AddLevel(std::make_unique<FiniteElementSpace>(fine_mesh, &rt_fec, attr_list));
    nd_h_u.AddLevel(std::make_unique<FiniteElementSpace>(mesh, &nd_fec));
    nd_h_u.AddLevel(std::make_unique<FiniteElementSpace>(fine_mesh, &nd_fec));
    CHECK(nd_h.GetFESpaceAtLevel(1).GlobalTrueVSize() > nd_fine.GlobalTrueVSize());
    const Operator &Ph_nd = nd_h.GetProlongationAtLevel(0);
    const Operator &Ph_rt = rt_h.GetProlongationAtLevel(0);
    CHECK(dynamic_cast<const mfem::HypreParMatrix *>(&Ph_nd));
    CheckCopies(nd_h.GetFESpaceAtLevel(0), nd_h.GetFESpaceAtLevel(1),
                nd_h_u.GetProlongationAtLevel(0), Ph_nd);
    const Operator &C_fine =
        rt_h.GetFESpaceAtLevel(1).GetDiscreteInterpolator(nd_h.GetFESpaceAtLevel(1));
    const Operator &C_c =
        rt_h.GetFESpaceAtLevel(0).GetDiscreteInterpolator(nd_h.GetFESpaceAtLevel(0));
    Vector x = Random(Ph_nd.Width()), px(Ph_nd.Height()), cpx(C_fine.Height()),
           cx(C_c.Height()), pcx(Ph_rt.Height());
    Ph_nd.Mult(x, px);
    C_fine.Mult(px, cpx);
    C_c.Mult(x, cx);
    Ph_rt.Mult(cx, pcx);
    pcx -= cpx;
    CHECK(MaxAbs(pcx) <= 1.0e-12 * MaxAbs(cpx));
  }
}

TEST_CASE("Broken space scalar curl", "[brokenspace][Serial][Parallel]")
{
  const auto comm = MPI_COMM_WORLD;
  const int order = GENERATE(1, 2);
  const auto type = GENERATE(mfem::Element::TRIANGLE, mfem::Element::QUADRILATERAL);
  const bool hanging = GENERATE(false, true);
  fem::DefaultIntegrationOrder::p_trial = order;
  auto smesh = MakeCrackedCubeMesh(4, type);
  smesh.EnsureNCMesh(true);
  auto part = Partition(smesh, Mpi::Size(comm));
  auto mesh = MakeParMesh(comm, smesh, part);
  const std::vector<int> attributes = {crack_attr};
  if (hanging)
  {
    RefineWhere(mesh, [](const mfem::Vector &c) { return c(0) > x_crack && c(1) < 0.5; });
    mesh.Update();
    REQUIRE(HasHangingEntities(mesh.Get()));
  }
  mfem::H1_FECollection h1_fec(order, 2);
  mfem::ND_FECollection nd_fec(order, 2);
  mfem::L2_FECollection l2_fec(order - 1, 2, mfem::BasisType::GaussLegendre,
                               mfem::FiniteElement::INTEGRAL);
  FiniteElementSpace h1(mesh, &h1_fec, attributes), nd(mesh, &nd_fec, attributes),
      l2(mesh, &l2_fec, attributes);
  const auto &gradient = nd.GetDiscreteInterpolator(h1);
  const auto &curl = l2.GetDiscreteInterpolator(nd);
  CHECK(h1.GlobalTrueVSize() > FiniteElementSpace(mesh, &h1_fec).GlobalTrueVSize());
  Vector x(h1.GetTrueVSize()), g(nd.GetTrueVSize()), cg(l2.GetTrueVSize());
  for (int i = 0; i < x.Size(); i++)
  {
    x(i) = std::sin(1.3 * (i + 1 + Mpi::Rank(comm)));
  }
  gradient.Mult(x, g);
  curl.Mult(g, cg);
  double error = cg.Normlinf(), scale = g.Normlinf();
  Mpi::GlobalMax(1, &error, comm);
  Mpi::GlobalMax(1, &scale, comm);
  CHECK(error < 1.0e-12 * scale);

  FiniteElementSpace nd_u(mesh, &nd_fec), l2_u(mesh, &l2_fec);
  Vector xu(nd_u.GetTrueVSize()), xb(nd.GetTrueVSize()), yu(l2_u.GetTrueVSize()),
      yb(l2.GetTrueVSize()), expected(l2.GetTrueVSize());
  for (int i = 0; i < xu.Size(); i++)
  {
    xu(i) = std::cos(0.7 * (i + 2 + Mpi::Rank(comm)));
  }
  nd.CopyFromConforming(xu, xb);
  curl.Mult(xb, yb);
  l2_u.GetDiscreteInterpolator(nd_u).Mult(xu, yu);
  l2.CopyFromConforming(yu, expected);
  scale = expected.Normlinf();
  yb -= expected;
  error = yb.Normlinf();
  Mpi::GlobalMax(1, &error, comm);
  Mpi::GlobalMax(1, &scale, comm);
  CHECK(scale > 1.0e-6);
  CHECK(error < 1.0e-12 * scale);
}

TEST_CASE("Broken hanging traces need more than row selection",
          "[brokenspace][Serial][Parallel]")
{
  const auto comm = MPI_COMM_WORLD;
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  const bool nonconforming = GENERATE(false, true);
  auto smesh = MakeCrackedCubeMesh(2, type);
  smesh.EnsureNCMesh(true);
  auto part = Partition(smesh, Mpi::Size(comm));
  auto mesh = MakeParMesh(comm, smesh, part);
  if (nonconforming)
  {
    RefineWhere(mesh, [](const mfem::Vector &c) { return c(0) > x_crack; });
    mesh.Update();
    REQUIRE(HasHangingEntities(mesh.Get()));
  }
  mfem::ND_FECollection nd_fec(2, 3);
  mfem::RT_FECollection rt_fec(1, 3);
  for (const mfem::FiniteElementCollection *fec :
       std::vector<const mfem::FiniteElementCollection *>{&nd_fec, &rt_fec})
  {
    FiniteElementSpace base(mesh, fec);
    FiniteElementSpace broken(mesh, fec, {crack_attr});
    const auto &P = *broken.Dof_TrueDof_Matrix();
    REQUIRE(P.GetGlobalNumCols() > base.GlobalTrueVSize());
    P.HostRead();
    mfem::SparseMatrix diag, offd;
    HYPRE_BigInt *cmap;
    P.GetDiag(diag);
    P.GetOffd(offd, cmap);

    // A selection left inverse can recover t only if a row of P is e_t^T.
    // Test the weaker necessary condition: a row with ANY nonzero multiple
    // of e_t^T, on ANY rank, even including rows not read by an element.
    // If this fails, changing ownership or allowing signed/scaled selection
    // cannot help. Include off-rank columns when counting row support.
    Vector singleton(P.Height()), available(P.Width());
    auto *mask = singleton.HostWrite();
    for (int l = 0; l < P.Height(); l++)
    {
      int nnz = 0;
      for (int k = diag.GetI()[l]; k < diag.GetI()[l + 1]; k++)
      {
        nnz += (std::abs(diag.GetData()[k]) > 1.0e-13);
      }
      for (int k = offd.GetI()[l]; k < offd.GetI()[l + 1]; k++)
      {
        nnz += (std::abs(offd.GetData()[k]) > 1.0e-13);
      }
      mask[l] = (nnz == 1);
    }
    P.AbsMultTranspose(1.0, singleton, 0.0, available);
    int counts[2] = {0, P.Width()};
    const auto *values = available.HostRead();
    for (int t = 0; t < P.Width(); t++)
    {
      counts[0] += (values[t] < 1.0e-13);
    }
    Mpi::GlobalSum(2, counts, comm);
    REQUIRE(counts[1] > 0);
    INFO(fec->Name() << ", type " << type << ", NC " << nonconforming << ": " << counts[0]
                     << " of " << counts[1]
                     << " true coefficients have no singleton row on any rank");
    if (nonconforming)
    {
      CHECK(counts[0] > 0);
    }
    else if (Mpi::Size(comm) == 1 || type == mfem::Element::HEXAHEDRON || fec == &rt_fec)
    {
      CHECK(counts[0] == 0);
    }
    else
    {
      // Parallel ND face transformations can already remove singleton rows
      // on a conforming sheet. Their count depends on the partition/orientation.
      CHECK(counts[0] <= counts[1]);
    }
  }
}

TEST_CASE("Broken space error estimators",
          "[brokenspace][errorestimator][Serial][Parallel]")
{
  const auto comm = MPI_COMM_WORLD;
  const int order = GENERATE_COPY(from_range(SerialOrParallel<int>(comm, {1, 2}, {2})));
  const auto type = GENERATE_COPY(from_range(SerialOrParallel<mfem::Element::Type>(
      comm, {mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON, mfem::Element::PYRAMID},
      {mfem::Element::TETRAHEDRON})));
  fem::DefaultIntegrationOrder::p_trial = order;

  config::MaterialData material;
  material.attributes = {1};
  material.epsilon_r.s = {2.0, 3.0, 4.0};
  material.mu_r.s = {1.0, 1.5, 2.0};
  config::PeriodicBoundaryData periodic;

  auto smesh = MakeCrackedCubeMesh(4, type);
  auto cut = CutMesh(smesh);
  auto part = Partition(smesh, Mpi::Size(comm));
  EstimatorSetup uncut_setup(comm, smesh, part, order, material, periodic);
  EstimatorSetup cut_setup(comm, cut, part, order, material, periodic);
  const std::vector<int> attr_list = {crack_attr}, no_attr_list = {};

  SECTION("Equivalence with a cut mesh")
  {
    // Broken recovery on the uncut mesh is the recovery on the cut mesh (on which the
    // interior boundary attribute is ignored, since it is not interior).
    const auto [grad_broken, curl_broken] =
        uncut_setup.Estimate(SmoothE, SmoothB, attr_list);
    const auto [grad_cut, curl_cut] = cut_setup.Estimate(SmoothE, SmoothB, attr_list);
    const auto [grad_cut_ref, curl_cut_ref] =
        cut_setup.Estimate(SmoothE, SmoothB, no_attr_list);
    CompareEstimates(comm, grad_broken, grad_cut, "gradient flux (broken vs. cut)");
    CompareEstimates(comm, curl_broken, curl_cut, "curl flux (broken vs. cut)");
    CompareEstimates(comm, grad_cut, grad_cut_ref, "gradient flux (cut)");
    CompareEstimates(comm, curl_cut, curl_cut_ref, "curl flux (cut)");

    // With multigrid requested, the broken spaces still use Jacobi preconditioning (without
    // a multigrid hierarchy), so the estimates are the same.
    const auto [grad_mg, curl_mg] = uncut_setup.Estimate(SmoothE, SmoothB, attr_list, true);
    CompareEstimates(comm, grad_mg, grad_broken, "gradient flux (multigrid requested)");
    CompareEstimates(comm, curl_mg, curl_broken, "curl flux (multigrid requested)");

    // Without the broken recovery the estimates differ.
    const auto [grad_cont, curl_cont] =
        uncut_setup.Estimate(SmoothE, SmoothB, no_attr_list);
    auto diff = grad_cont;
    diff -= grad_broken;
    CHECK(GlobalMax(comm, diff) > 1.0e-3 * GlobalMax(comm, grad_cont));
    diff = curl_cont;
    diff -= curl_broken;
    CHECK(GlobalMax(comm, diff) > 1.0e-3 * GlobalMax(comm, curl_cont));
  }

  SECTION("Exact recovery of a jump")
  {
    // A piecewise constant flux with a jump across the interior boundary is recovered
    // exactly by the broken spaces, but not by the continuous ones. (The recovery on the
    // mesh with pyramids is not exact for piecewise constant fields at higher orders, also
    // without interior boundaries.)
    if (type == mfem::Element::PYRAMID && order > 1)
    {
      return;
    }
    const auto [grad_broken, curl_broken] = uncut_setup.Estimate(JumpE, JumpB, attr_list);
    const auto [grad_cont, curl_cont] = uncut_setup.Estimate(JumpE, JumpB, no_attr_list);
    const double grad_max = GlobalMax(comm, grad_cont),
                 curl_max = GlobalMax(comm, curl_cont);
    CHECK(grad_max > 1.0e-2);
    CHECK(curl_max > 1.0e-2);
    CHECK(GlobalMax(comm, grad_broken) < 1.0e-6 * grad_max);
    CHECK(GlobalMax(comm, curl_broken) < 1.0e-6 * curl_max);
  }
}

namespace
{

// Element-wise values matched between two meshes of the same geometry and partitioning by
// sorting the local elements by their centers.
std::vector<double> SortByCenter(Mesh &mesh, const Vector &v)
{
  auto &pmesh = mesh.Get();
  std::vector<std::pair<std::array<double, 3>, double>> entries;
  mfem::Vector c(pmesh.SpaceDimension());
  for (int e = 0; e < pmesh.GetNE(); e++)
  {
    pmesh.GetElementCenter(e, c);
    std::array<double, 3> key = {0.0, 0.0, 0.0};
    for (int d = 0; d < c.Size(); d++)
    {
      key[d] = std::round(c(d) * 1.0e8) * 1.0e-8;
    }
    entries.emplace_back(key, v(e));
  }
  std::sort(entries.begin(), entries.end(),
            [](const auto &a, const auto &b) { return a.first < b.first; });
  std::vector<double> sorted;
  for (const auto &entry : entries)
  {
    sorted.push_back(entry.second);
  }
  return sorted;
}

double GlobalSumSquares(MPI_Comm comm, const Vector &v)
{
  double sum = v * v;
  Mpi::GlobalSum(1, &sum, comm);
  return sum;
}

// Complete a mesh modification: update the mesh and the finite element spaces.
void UpdateSetup(Mesh &mesh, FiniteElementSpaceHierarchy &nd_fespaces,
                 FiniteElementSpaceHierarchy &rt_fespaces)
{
  mesh.Update();
  for (auto *fespaces : {&nd_fespaces, &rt_fespaces})
  {
    fespaces->GetFinestFESpace().Get().Update(false);
    fespaces->GetFinestFESpace().Update();
  }
}

}  // namespace

TEST_CASE("Broken space error estimators (nonconforming)",
          "[brokenspace][errorestimator][Serial][Parallel]")
{
  const auto comm = MPI_COMM_WORLD;
  const int order = GENERATE_COPY(from_range(SerialOrParallel<int>(comm, {1, 2}, {2})));
  const auto type = GENERATE_COPY(from_range(SerialOrParallel<mfem::Element::Type>(
      comm, {mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON},
      {mfem::Element::HEXAHEDRON})));
  fem::DefaultIntegrationOrder::p_trial = order;

  config::MaterialData material;
  material.attributes = {1};
  material.epsilon_r.s = {2.0, 3.0, 4.0};
  material.mu_r.s = {1.0, 1.5, 2.0};
  config::PeriodicBoundaryData periodic;

  auto smesh = MakeCrackedCubeMesh(4, type);
  auto cut = CutMesh(smesh);
  smesh.EnsureNCMesh(true);
  cut.EnsureNCMesh(true);
  auto part = Partition(smesh, Mpi::Size(comm));
  auto cut_part = Partition(cut, Mpi::Size(comm));
  EstimatorSetup uncut_setup(comm, smesh, part, order, material, periodic);
  EstimatorSetup cut_setup(comm, cut, cut_part, order, material, periodic);
  const std::vector<int> attr_list = {crack_attr}, no_attr_list = {};

  // Sides are computed before the refinement (as by an estimator in an adaptive loop).

  SECTION("Equivalence with a cut mesh")
  {
    // Refinement which is the same on both sides of the interior boundary, and ends on it
    // (hanging entities on the interior boundary and away from it).
    for (auto *setup : {&uncut_setup, &cut_setup})
    {
      RefineWhere(setup->mesh, [](const mfem::Vector &c)
                  { return std::abs(c(0) - x_crack) < 0.25 && c(1) < 0.5; });
      UpdateSetup(setup->mesh, setup->nd_fespaces, setup->rt_fespaces);
    }
    REQUIRE(HasHangingEntities(uncut_setup.mesh.Get()));
    const auto [grad_broken, curl_broken] =
        uncut_setup.Estimate(SmoothE, SmoothB, attr_list);
    const auto [grad_cut, curl_cut] = cut_setup.Estimate(SmoothE, SmoothB, no_attr_list);
    auto Compare = [&](const Vector &a, const Vector &b, const char *what)
    {
      const auto sa = SortByCenter(uncut_setup.mesh, a);
      const auto sb = SortByCenter(cut_setup.mesh, b);
      CHECK(sa.size() == sb.size());  // Not REQUIRE: the reductions below are collective
      double max_diff = 0.0, max_ref = 0.0;
      for (std::size_t i = 0; i < std::min(sa.size(), sb.size()); i++)
      {
        max_diff = std::max(max_diff, std::abs(sa[i] - sb[i]));
        max_ref = std::max(max_ref, std::abs(sb[i]));
      }
      Mpi::GlobalMax(1, &max_diff, comm);
      Mpi::GlobalMax(1, &max_ref, comm);
      INFO(what << ": max. difference " << max_diff << " (max. estimate " << max_ref
                << ")");
      CHECK(max_diff <= 1.0e-8 * max_ref + 1.0e-14);
    };
    Compare(grad_broken, grad_cut, "gradient flux (broken vs. cut)");
    Compare(curl_broken, curl_cut, "curl flux (broken vs. cut)");

    // With multigrid requested, the broken spaces still use Jacobi preconditioning (without
    // a multigrid hierarchy), so the estimates are the same.
    const auto [grad_mg, curl_mg] = uncut_setup.Estimate(SmoothE, SmoothB, attr_list, true);
    CompareEstimates(comm, grad_mg, grad_broken, "gradient flux (multigrid requested)");
    CompareEstimates(comm, curl_mg, curl_broken, "curl flux (multigrid requested)");
  }

  SECTION("Exact recovery of a jump with hanging entities across the interior boundary")
  {
    // Refinement of one side only, then of part of the other side.
    RefineWhere(uncut_setup.mesh, [](const mfem::Vector &c)
                { return c(0) < x_crack && c(0) > 0.25 && c(2) < 0.75; });
    UpdateSetup(uncut_setup.mesh, uncut_setup.nd_fespaces, uncut_setup.rt_fespaces);
    RefineWhere(uncut_setup.mesh, [](const mfem::Vector &c)
                { return c(0) > x_crack && c(0) < 0.75 && c(1) < 0.4; });
    UpdateSetup(uncut_setup.mesh, uncut_setup.nd_fespaces, uncut_setup.rt_fespaces);
    REQUIRE(HasHangingEntities(uncut_setup.mesh.Get()));
    const auto [grad_broken, curl_broken] = uncut_setup.Estimate(JumpE, JumpB, attr_list);
    const auto [grad_cont, curl_cont] = uncut_setup.Estimate(JumpE, JumpB, no_attr_list);
    const double grad_max = GlobalMax(comm, grad_cont),
                 curl_max = GlobalMax(comm, curl_cont);
    CHECK(grad_max > 1.0e-4);
    CHECK(curl_max > 1.0e-4);
    CHECK(GlobalMax(comm, grad_broken) < 1.0e-6 * grad_max);
    CHECK(GlobalMax(comm, curl_broken) < 1.0e-6 * curl_max);
  }

  SECTION("Interior boundary sides through rebalancing")
  {
    // Two identical setups, one of which is rebalanced after the refinement.
    if (Mpi::Size(comm) > 1)
    {
      EstimatorSetup other_setup(comm, smesh, part, order, material, periodic);
      for (auto *setup : {&uncut_setup, &other_setup})
      {
        RefineWhere(setup->mesh, [](const mfem::Vector &c)
                    { return c(0) < x_crack && c(0) > 0.25 && c(2) < 0.75; });
        UpdateSetup(setup->mesh, setup->nd_fespaces, setup->rt_fespaces);
        RefineWhere(setup->mesh, [](const mfem::Vector &c)
                    { return c(0) > x_crack && c(0) < 0.75 && c(1) < 0.4; });
        if (setup == &uncut_setup)
        {
          setup->mesh.Get().Rebalance();
        }
        UpdateSetup(setup->mesh, setup->nd_fespaces, setup->rt_fespaces);
      }
      const auto [grad_reb, curl_reb] = uncut_setup.Estimate(SmoothE, SmoothB, attr_list);
      const auto [grad_ref, curl_ref] = other_setup.Estimate(SmoothE, SmoothB, attr_list);
      const auto [grad_cont, curl_cont] =
          other_setup.Estimate(SmoothE, SmoothB, no_attr_list);
      CHECK_THAT(GlobalSumSquares(comm, grad_reb),
                 WithinRel(GlobalSumSquares(comm, grad_ref), 1.0e-8));
      CHECK_THAT(GlobalSumSquares(comm, curl_reb),
                 WithinRel(GlobalSumSquares(comm, curl_ref), 1.0e-8));
      CHECK(std::abs(GlobalSumSquares(comm, grad_cont) - GlobalSumSquares(comm, grad_ref)) >
            1.0e-3 * GlobalSumSquares(comm, grad_ref));
    }
  }
}

TEST_CASE("Broken space error estimators (nonconforming input mesh)",
          "[brokenspace][errorestimator][Serial][Parallel]")
{
  // A nonconforming input mesh, refined before the simulation, with pyramids or with
  // anisotropic refinements, whose interior boundary sides are inherited from those of its
  // root elements. With the same refinement on both sides of the interior boundary, the
  // estimates with the broken spaces are those of the cut mesh refined in the same way.
  const auto comm = MPI_COMM_WORLD;
  const auto type = GENERATE_COPY(from_range(SerialOrParallel<mfem::Element::Type>(
      comm, {mfem::Element::PYRAMID, mfem::Element::HEXAHEDRON, mfem::Element::TETRAHEDRON},
      {mfem::Element::PYRAMID, mfem::Element::HEXAHEDRON})));
  CAPTURE(type);
  constexpr int order = 1;
  fem::DefaultIntegrationOrder::p_trial = order;
  config::MaterialData material;
  material.attributes = {1};
  material.epsilon_r.s = {2.0, 3.0, 4.0};
  material.mu_r.s = {1.0, 1.5, 2.0};
  config::PeriodicBoundaryData periodic;

  // Two levels of refinement of the elements touching the interior boundary (for y < 0.5),
  // anisotropic (normal to it) for the mesh of hexahedra, and of the tetrahedral children
  // of the pyramids. The refinement is the same on both sides of the interior boundary (the
  // estimates differ where only one side is refined). MFEM does not refine meshes with
  // pyramids and nodes, so the meshes are refined before creating the setups.
  auto smesh = MakeCrackedCubeMesh(4, type);
  auto cut = CutMesh(smesh);
  for (auto *m : {&smesh, &cut})
  {
    m->EnsureNCMesh(true);
    for (int level = 0; level < 2; level++)
    {
      mfem::Array<mfem::Refinement> refs;
      mfem::Array<int> v;
      mfem::Vector c;
      for (int e = 0; e < m->GetNE(); e++)
      {
        m->GetElementCenter(e, c);
        m->GetElementVertices(e, v);
        if (c(1) < 0.5 &&
            std::any_of(v.begin(), v.end(), [&](int i)
                        { return std::abs(m->GetVertex(i)[0] - x_crack) < 1.0e-12; }))
        {
          refs.Append((type == mfem::Element::HEXAHEDRON)
                          ? mfem::Refinement(e, mfem::Refinement::X)
                          : mfem::Refinement(e));
        }
      }
      m->GeneralRefinement(refs, 1, 0);
    }
  }
  auto part = Partition(smesh, Mpi::Size(comm));
  auto cut_part = Partition(cut, Mpi::Size(comm));
  EstimatorSetup uncut_setup(comm, smesh, part, order, material, periodic);
  EstimatorSetup cut_setup(comm, cut, cut_part, order, material, periodic);
  REQUIRE(HasHangingEntities(uncut_setup.mesh.Get()));
  const std::vector<int> attr_list = {crack_attr}, no_attr_list = {};
  mfem::ND_FECollection nd(1, 3);
  FiniteElementSpace broken(uncut_setup.mesh, &nd, attr_list),
      continuous(uncut_setup.mesh, &nd);
  REQUIRE(broken.GlobalTrueVSize() > continuous.GlobalTrueVSize());

  const auto [grad_broken, curl_broken] = uncut_setup.Estimate(SmoothE, SmoothB, attr_list);
  const auto [grad_cut, curl_cut] = cut_setup.Estimate(SmoothE, SmoothB, no_attr_list);
  auto Compare = [&](const Vector &a, const Vector &b, const char *what)
  {
    const auto sa = SortByCenter(uncut_setup.mesh, a);
    const auto sb = SortByCenter(cut_setup.mesh, b);
    CHECK(sa.size() == sb.size());  // Not REQUIRE: the reductions below are collective
    double max_diff = 0.0, max_ref = 0.0;
    for (std::size_t i = 0; i < std::min(sa.size(), sb.size()); i++)
    {
      max_diff = std::max(max_diff, std::abs(sa[i] - sb[i]));
      max_ref = std::max(max_ref, std::abs(sb[i]));
    }
    Mpi::GlobalMax(1, &max_diff, comm);
    Mpi::GlobalMax(1, &max_ref, comm);
    INFO(what << ": max. difference " << max_diff << " (max. estimate " << max_ref << ")");
    CHECK(max_diff <= 1.0e-8 * max_ref + 1.0e-14);
  };
  Compare(grad_broken, grad_cut, "gradient flux (broken vs. cut)");
  Compare(curl_broken, curl_cut, "curl flux (broken vs. cut)");
}

TEST_CASE("Broken space error estimators (saved nonconforming mesh)",
          "[brokenspace][errorestimator][Serial][Parallel]")
{
  // A nonconforming mesh refined in a previous simulation, saved and loaded as for a
  // restart: the interior boundary sides are inherited from the root elements of its
  // refinement hierarchy.
  const auto comm = MPI_COMM_WORLD;
  const int order = GENERATE_COPY(from_range(SerialOrParallel<int>(comm, {1, 2}, {1})));
  const auto type = GENERATE_COPY(from_range(SerialOrParallel<mfem::Element::Type>(
      comm, {mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON},
      {mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON})));
  fem::DefaultIntegrationOrder::p_trial = order;
  config::MaterialData material;
  material.attributes = {1};
  material.epsilon_r.s = {2.0, 3.0, 4.0};
  material.mu_r.s = {1.0, 1.5, 2.0};
  config::PeriodicBoundaryData periodic;

  auto smesh = MakeCrackedCubeMesh(4, type);
  smesh.EnsureNCMesh(true);
  auto part = Partition(smesh, Mpi::Size(comm));
  EstimatorSetup setup(comm, smesh, part, order, material, periodic);
  const std::vector<int> attr_list = {crack_attr}, no_attr_list = {};

  // Two levels of refinement: of one side, then across the interior boundary.
  RefineWhere(setup.mesh, [](const mfem::Vector &c)
              { return c(0) < x_crack && c(0) > 0.25 && c(2) < 0.75; });
  UpdateSetup(setup.mesh, setup.nd_fespaces, setup.rt_fespaces);
  RefineWhere(setup.mesh, [](const mfem::Vector &c)
              { return std::abs(c(0) - x_crack) < 0.2 && c(1) < 0.4; });
  UpdateSetup(setup.mesh, setup.nd_fespaces, setup.rt_fespaces);
  REQUIRE(HasHangingEntities(setup.mesh.Get()));

  // Save the mesh as for config["Model"]["Refinement"]["SaveAdaptMesh"], and load it on all
  // processes.
  std::string saved;
  {
    mfem::ParMesh gmesh(setup.mesh.Get());
    mfem::Array<int> serial_partition(gmesh.GetNE());
    serial_partition = 0;
    gmesh.Rebalance(serial_partition);
    if (Mpi::Root(comm))
    {
      std::ostringstream fo;
      fo.precision(17);
      gmesh.Mesh::Print(fo);
      saved = fo.str();
    }
    int len = static_cast<int>(saved.size());
    Mpi::Broadcast(1, &len, 0, comm);
    saved.resize(len);
    Mpi::Broadcast(len, saved.data(), 0, comm);
  }
  std::istringstream fi(saved);
  mfem::Mesh loaded(fi, 0, 1, false);
  REQUIRE(loaded.Nonconforming());
  auto loaded_part = Partition(loaded, Mpi::Size(comm));
  EstimatorSetup loaded_setup(comm, loaded, loaded_part, order, material, periodic);
  REQUIRE(HasHangingEntities(loaded_setup.mesh.Get()));

  // The partitioning keeps the elements of the loaded mesh on the processes of the refined
  // one, so that the estimates can be compared locally.
  const auto [grad_ref, curl_ref] = setup.Estimate(SmoothE, SmoothB, attr_list);
  const auto [grad_loaded, curl_loaded] =
      loaded_setup.Estimate(SmoothE, SmoothB, attr_list);
  const auto [grad_cont, curl_cont] = loaded_setup.Estimate(SmoothE, SmoothB, no_attr_list);
  auto Compare = [&](const Vector &a, const Vector &b, const char *what)
  {
    const auto sa = SortByCenter(loaded_setup.mesh, a);
    const auto sb = SortByCenter(setup.mesh, b);
    CHECK(sa.size() == sb.size());  // Not REQUIRE: the reductions below are collective
    double max_diff = 0.0, max_ref = 0.0;
    for (std::size_t i = 0; i < std::min(sa.size(), sb.size()); i++)
    {
      max_diff = std::max(max_diff, std::abs(sa[i] - sb[i]));
      max_ref = std::max(max_ref, std::abs(sb[i]));
    }
    Mpi::GlobalMax(1, &max_diff, comm);
    Mpi::GlobalMax(1, &max_ref, comm);
    INFO(what << ": max. difference " << max_diff << " (max. estimate " << max_ref << ")");
    CHECK(max_diff <= 1.0e-8 * max_ref + 1.0e-14);
  };
  Compare(grad_loaded, grad_ref, "gradient flux (loaded vs. refined)");
  Compare(curl_loaded, curl_ref, "curl flux (loaded vs. refined)");
  CHECK(std::abs(GlobalSumSquares(comm, grad_cont) - GlobalSumSquares(comm, grad_ref)) >
        1.0e-3 * GlobalSumSquares(comm, grad_ref));
}

TEST_CASE("Broken space error estimators (two interior boundaries)",
          "[brokenspace][errorestimator][Serial][Parallel]")
{
  // Two parallel interior boundaries with a single layer of elements between them. With a
  // nonconforming refinement of part of the layer, the DOFs of hanging entities are
  // constrained by those of both interior boundaries, so the elements of the layer have to
  // read either the original DOFs of both or the copies of both. The recovery of a
  // piecewise constant field with jumps across both boundaries is then exact, independently
  // of the numbering of the elements and of the partitioning.
  const auto comm = MPI_COMM_WORLD;
  const int order = GENERATE_COPY(from_range(SerialOrParallel<int>(comm, {1, 2}, {1})));
  const auto type = GENERATE_COPY(from_range(SerialOrParallel<mfem::Element::Type>(
      comm, {mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON},
      {mfem::Element::HEXAHEDRON})));
  const bool nonconforming =
      GENERATE_COPY(from_range(SerialOrParallel<int>(comm, {0, 1}, {1})));
  const bool layer_first = GENERATE(false, true);
  const bool swap_ranks =
      GENERATE_COPY(from_range(SerialOrParallel<int>(comm, {0}, {0, 1})));
  fem::DefaultIntegrationOrder::p_trial = order;
  config::MaterialData material;
  material.attributes = {1};
  material.epsilon_r.s = {2.0, 3.0, 4.0};
  material.mu_r.s = {1.0, 1.5, 2.0};
  config::PeriodicBoundaryData periodic;

  auto InLayer = [](const mfem::Vector &c) { return c(0) > x_sheet && c(0) < x_crack; };
  auto smesh = MakeCrackedCubeMesh(4, type, 1.0, {x_sheet, x_crack});
  if (layer_first)
  {
    smesh = ReorderElements(smesh, InLayer);
  }
  if (nonconforming)
  {
    smesh.EnsureNCMesh(true);
  }
  auto part = Partition(smesh, Mpi::Size(comm));
  if (swap_ranks)
  {
    for (auto &p : part)
    {
      p = Mpi::Size(comm) - 1 - p;
    }
  }
  EstimatorSetup setup(comm, smesh, part, order, material, periodic);
  const std::vector<int> attr_list = {crack_attr}, no_attr_list = {};
  if (nonconforming)
  {
    RefineWhere(setup.mesh,
                [&](const mfem::Vector &c) { return InLayer(c) && c(1) < 0.5; });
    UpdateSetup(setup.mesh, setup.nd_fespaces, setup.rt_fespaces);
    REQUIRE(HasHangingEntities(setup.mesh.Get()));
  }
  const auto [grad_broken, curl_broken] = setup.Estimate(JumpE2, JumpB2, attr_list);
  const auto [grad_cont, curl_cont] = setup.Estimate(JumpE2, JumpB2, no_attr_list);
  const double grad_max = GlobalMax(comm, grad_cont), curl_max = GlobalMax(comm, curl_cont);
  CAPTURE(order, type, nonconforming, layer_first, swap_ranks);
  CHECK(grad_max > 1.0e-4);
  CHECK(curl_max > 1.0e-4);
  CHECK(GlobalMax(comm, grad_broken) < 1.0e-6 * grad_max);
  CHECK(GlobalMax(comm, curl_broken) < 1.0e-6 * curl_max);
}

TEST_CASE("Broken space error estimators (junctions)",
          "[brokenspace][errorestimator][Serial][Parallel]")
{
  // A tunnel through a ground plane: its legs meet the ground plane at junctions of three
  // regions, each of which reads its own version of the DOFs of the junction edges. The
  // recovery of piecewise constant fields with different values in the three regions is
  // exact, also with a nonconforming refinement next to a junction on one side only (with
  // hanging entities constrained by the DOFs of the junction edges).
  const auto comm = MPI_COMM_WORLD;
  const int order = GENERATE_COPY(from_range(SerialOrParallel<int>(comm, {1, 2}, {2})));
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  const bool nonconforming = GENERATE(false, true);
  fem::DefaultIntegrationOrder::p_trial = order;
  config::MaterialData material;
  material.attributes = {1};
  material.epsilon_r.s = {2.0, 3.0, 4.0};
  material.mu_r.s = {1.0, 1.5, 2.0};
  config::PeriodicBoundaryData periodic;

  auto smesh = MakeAirBridgeMesh(4, type, 0.0, 1.0, 0.75);
  if (nonconforming)
  {
    smesh.EnsureNCMesh(true);
  }
  auto part = Partition(smesh, Mpi::Size(comm));
  EstimatorSetup setup(comm, smesh, part, order, material, periodic);
  const std::vector<int> attr_list = {crack_attr}, no_attr_list = {};
  if (nonconforming)
  {
    // The elements in the tunnel next to the junction x = 0.25, z = 0.25.
    RefineWhere(setup.mesh, [](const mfem::Vector &c)
                { return c(0) > 0.25 && c(0) < 0.5 && c(2) > 0.25 && c(2) < 0.5; });
    UpdateSetup(setup.mesh, setup.nd_fespaces, setup.rt_fespaces);
    REQUIRE(HasHangingEntities(setup.mesh.Get()));
  }
  const auto [grad_broken, curl_broken] =
      setup.Estimate(JumpETunnel, JumpBTunnel, attr_list);
  const auto [grad_cont, curl_cont] =
      setup.Estimate(JumpETunnel, JumpBTunnel, no_attr_list);
  const double grad_max = GlobalMax(comm, grad_cont), curl_max = GlobalMax(comm, curl_cont);
  CAPTURE(order, type, nonconforming);
  CHECK(grad_max > 1.0e-4);
  CHECK(curl_max > 1.0e-4);
  CHECK(GlobalMax(comm, grad_broken) < 1.0e-6 * grad_max);
  CHECK(GlobalMax(comm, curl_broken) < 1.0e-6 * curl_max);
}

TEST_CASE("Broken space error estimators (partial interior boundary)",
          "[brokenspace][errorestimator][Serial][Parallel]")
{
  // An interior boundary with a free edge (at y = 0.5), where the recovery remains
  // continuous, and a partitioning where some processes do not touch the interior boundary.
  const auto comm = MPI_COMM_WORLD;
  const int order = GENERATE_COPY(from_range(SerialOrParallel<int>(comm, {1, 2}, {1})));
  const auto type = GENERATE_COPY(from_range(SerialOrParallel<mfem::Element::Type>(
      comm, {mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON},
      {mfem::Element::TETRAHEDRON})));
  const bool nonconforming = GENERATE(false, true);
  fem::DefaultIntegrationOrder::p_trial = order;
  config::MaterialData material;
  material.attributes = {1};
  material.epsilon_r.s = {2.0, 3.0, 4.0};
  material.mu_r.s = {1.0, 1.5, 2.0};
  config::PeriodicBoundaryData periodic;
  constexpr double y_max = 0.5;

  auto smesh = MakeCrackedCubeMesh(4, type, y_max);
  auto cut = CutMesh(smesh, y_max);
  if (nonconforming)
  {
    smesh.EnsureNCMesh(true);
    cut.EnsureNCMesh(true);
  }
  auto PartitionY = [](const mfem::Mesh &mesh, int num_procs)
  {
    std::vector<int> part(mesh.GetNE());
    mfem::Vector c(3);
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      const_cast<mfem::Mesh &>(mesh).GetElementCenter(e, c);
      part[e] = std::min(num_procs - 1, static_cast<int>(c(1) * num_procs));
    }
    return part;
  };
  auto part = PartitionY(smesh, Mpi::Size(comm));
  auto cut_part = PartitionY(cut, Mpi::Size(comm));
  EstimatorSetup uncut_setup(comm, smesh, part, order, material, periodic);
  EstimatorSetup cut_setup(comm, cut, cut_part, order, material, periodic);
  const std::vector<int> attr_list = {crack_attr}, no_attr_list = {};
  if (nonconforming)
  {
    auto pred = [](const mfem::Vector &c) { return c(0) < 0.75 && c(1) < 0.75; };
    RefineWhere(uncut_setup.mesh, pred);
    UpdateSetup(uncut_setup.mesh, uncut_setup.nd_fespaces, uncut_setup.rt_fespaces);
    RefineWhere(cut_setup.mesh, pred);
    UpdateSetup(cut_setup.mesh, cut_setup.nd_fespaces, cut_setup.rt_fespaces);
  }
  const auto [grad_broken, curl_broken] = uncut_setup.Estimate(SmoothE, SmoothB, attr_list);
  const auto [grad_cut, curl_cut] = cut_setup.Estimate(SmoothE, SmoothB, no_attr_list);
  auto Compare = [&](const Vector &a, const Vector &b, const char *what)
  {
    const auto sa = SortByCenter(uncut_setup.mesh, a);
    const auto sb = SortByCenter(cut_setup.mesh, b);
    CHECK(sa.size() == sb.size());
    double max_diff = 0.0, max_ref = 0.0;
    for (std::size_t i = 0; i < std::min(sa.size(), sb.size()); i++)
    {
      max_diff = std::max(max_diff, std::abs(sa[i] - sb[i]));
      max_ref = std::max(max_ref, std::abs(sb[i]));
    }
    Mpi::GlobalMax(1, &max_diff, comm);
    Mpi::GlobalMax(1, &max_ref, comm);
    INFO(what << ": max. difference " << max_diff << " (max. estimate " << max_ref << ")");
    CHECK(max_diff <= 1.0e-8 * max_ref + 1.0e-14);
  };
  Compare(grad_broken, grad_cut, "gradient flux (broken vs. cut)");
  Compare(curl_broken, curl_cut, "curl flux (broken vs. cut)");
}

TEST_CASE("Broken space error estimators (2D)",
          "[brokenspace][errorestimator][Serial][Parallel]")
{
  const auto comm = MPI_COMM_WORLD;
  const int order = GENERATE_COPY(from_range(SerialOrParallel<int>(comm, {1, 2}, {2})));
  const auto type = GENERATE_COPY(from_range(SerialOrParallel<mfem::Element::Type>(
      comm, {mfem::Element::TRIANGLE, mfem::Element::QUADRILATERAL},
      {mfem::Element::TRIANGLE})));
  fem::DefaultIntegrationOrder::p_trial = order;

  config::MaterialData material;
  material.attributes = {1};
  material.epsilon_r.s = {2.0, 3.0, 4.0};
  material.mu_r.s = {1.0, 1.5, 2.0};
  config::PeriodicBoundaryData periodic;

  auto smesh = MakeCrackedCubeMesh(6, type);
  auto cut = CutMesh(smesh);
  auto part = Partition(smesh, Mpi::Size(comm));
  EstimatorSetup2D uncut_setup(comm, smesh, part, order, material, periodic);
  EstimatorSetup2D cut_setup(comm, cut, part, order, material, periodic);
  const std::vector<int> attr_list = {crack_attr}, no_attr_list = {};

  SECTION("Equivalence with a cut mesh")
  {
    const auto [grad_broken, curl_broken] =
        uncut_setup.Estimate(SmoothE2D, SmoothB2D, attr_list);
    const auto [grad_cut, curl_cut] =
        cut_setup.Estimate(SmoothE2D, SmoothB2D, no_attr_list);
    CompareEstimates(comm, grad_broken, grad_cut, "gradient flux (broken vs. cut)");
    CompareEstimates(comm, curl_broken, curl_cut, "curl flux (broken vs. cut)");
    const auto [grad_cont, curl_cont] =
        uncut_setup.Estimate(SmoothE2D, SmoothB2D, no_attr_list);
    auto diff = grad_cont;
    diff -= grad_broken;
    CHECK(GlobalMax(comm, diff) > 1.0e-3 * GlobalMax(comm, grad_cont));
    diff = curl_cont;
    diff -= curl_broken;
    CHECK(GlobalMax(comm, diff) > 1.0e-3 * GlobalMax(comm, curl_cont));
  }

  SECTION("Exact recovery of a jump")
  {
    const auto [grad_broken, curl_broken] =
        uncut_setup.Estimate(JumpE2D, JumpB2D, attr_list);
    const auto [grad_cont, curl_cont] =
        uncut_setup.Estimate(JumpE2D, JumpB2D, no_attr_list);
    const double grad_max = GlobalMax(comm, grad_cont),
                 curl_max = GlobalMax(comm, curl_cont);
    CHECK(grad_max > 1.0e-4);
    CHECK(curl_max > 1.0e-4);
    CHECK(GlobalMax(comm, grad_broken) < 1.0e-6 * grad_max);
    CHECK(GlobalMax(comm, curl_broken) < 1.0e-6 * curl_max);
  }
}

TEST_CASE("Broken space error estimators (2D nonconforming mesh without history)",
          "[brokenspace][errorestimator][Serial][Parallel]")
{
  // A mesh refined before its interior boundary sides are computed, as one loaded from a
  // previous simulation: the sides are inherited from the root elements. The
  // refinement is the same on both sides of the interior boundary, for the equivalence with
  // a cut mesh (the limits in y are on element boundaries, so that the triangles on both
  // sides of an interior boundary edge are refined alike).
  const auto comm = MPI_COMM_WORLD;
  const int order = GENERATE_COPY(from_range(SerialOrParallel<int>(comm, {1, 2}, {1})));
  const auto type = GENERATE_COPY(from_range(SerialOrParallel<mfem::Element::Type>(
      comm, {mfem::Element::TRIANGLE, mfem::Element::QUADRILATERAL},
      {mfem::Element::QUADRILATERAL})));
  fem::DefaultIntegrationOrder::p_trial = order;
  config::MaterialData material;
  material.attributes = {1};
  material.epsilon_r.s = {2.0, 3.0, 4.0};
  material.mu_r.s = {1.0, 1.5, 2.0};
  config::PeriodicBoundaryData periodic;

  auto smesh = MakeCrackedCubeMesh(6, type);
  auto cut = CutMesh(smesh);
  for (auto *m : {&smesh, &cut})
  {
    m->EnsureNCMesh(true);
    for (const auto &[dx, y_max] : {std::pair{0.3, 0.5}, std::pair{0.15, 0.25}})
    {
      mfem::Array<int> marked;
      mfem::Vector c(2);
      for (int e = 0; e < m->GetNE(); e++)
      {
        m->GetElementCenter(e, c);
        if (std::abs(c(0) - x_crack) < dx && c(1) < y_max)
        {
          marked.Append(e);
        }
      }
      m->GeneralRefinement(marked, 1, 0);
    }
  }
  auto part = Partition(smesh, Mpi::Size(comm));
  auto cut_part = Partition(cut, Mpi::Size(comm));
  EstimatorSetup2D uncut_setup(comm, smesh, part, order, material, periodic);
  EstimatorSetup2D cut_setup(comm, cut, cut_part, order, material, periodic);
  REQUIRE(HasHangingEntities(uncut_setup.mesh.Get()));
  const std::vector<int> attr_list = {crack_attr}, no_attr_list = {};
  const auto [grad_broken, curl_broken] =
      uncut_setup.Estimate(SmoothE2D, SmoothB2D, attr_list);
  const auto [grad_cut, curl_cut] = cut_setup.Estimate(SmoothE2D, SmoothB2D, no_attr_list);
  const auto [grad_cont, curl_cont] =
      uncut_setup.Estimate(SmoothE2D, SmoothB2D, no_attr_list);
  auto Compare = [&](const Vector &a, const Vector &b, const char *what)
  {
    const auto sa = SortByCenter(uncut_setup.mesh, a);
    const auto sb = SortByCenter(cut_setup.mesh, b);
    CHECK(sa.size() == sb.size());  // Not REQUIRE: the reductions below are collective
    double max_diff = 0.0, max_ref = 0.0;
    for (std::size_t i = 0; i < std::min(sa.size(), sb.size()); i++)
    {
      max_diff = std::max(max_diff, std::abs(sa[i] - sb[i]));
      max_ref = std::max(max_ref, std::abs(sb[i]));
    }
    Mpi::GlobalMax(1, &max_diff, comm);
    Mpi::GlobalMax(1, &max_ref, comm);
    INFO(what << ": max. difference " << max_diff << " (max. estimate " << max_ref << ")");
    CHECK(max_diff <= 1.0e-8 * max_ref + 1.0e-14);
  };
  Compare(grad_broken, grad_cut, "gradient flux (broken vs. cut)");
  Compare(curl_broken, curl_cut, "curl flux (broken vs. cut)");
  CHECK(std::abs(GlobalSumSquares(comm, grad_cont) - GlobalSumSquares(comm, grad_broken)) >
        1.0e-3 * GlobalSumSquares(comm, grad_broken));
}

TEST_CASE("Interface dielectric energy on interior boundaries separating fields",
          "[brokenspace][surfacefunctional][Serial][Parallel]")
{
  // On an interior boundary which separates the fields on its two sides (such as a thin
  // metal sheet), the interface energy sums the energies of both sides, as for a mesh cut
  // along the boundary.
  const auto comm = MPI_COMM_WORLD;
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  const int order = GENERATE(1, 2);
  const bool dielectric = GENERATE(false, true);
  fem::DefaultIntegrationOrder::p_trial = order;
  config::MaterialData material;
  material.attributes = {1};
  if (dielectric)
  {
    material.epsilon_r.s = {11.7, 11.7, 11.7};
  }
  config::PeriodicBoundaryData periodic;

  auto smesh = MakeCrackedCubeMesh(4, type);
  auto cut = CutMesh(smesh);
  auto part = Partition(smesh, Mpi::Size(comm));
  auto cut_part = Partition(cut, Mpi::Size(comm));
  auto mesh = MakeParMesh(comm, smesh, part);
  auto cut_mesh = MakeParMesh(comm, cut, cut_part);
  MaterialOperator mat_op({material}, periodic, ProblemType::EIGENMODE, mesh);
  MaterialOperator cut_mat_op({material}, periodic, ProblemType::EIGENMODE, cut_mesh);
  mfem::ND_FECollection nd_fec(order, 3);
  FiniteElementSpace nd_fespace(mesh, &nd_fec), cut_nd_fespace(cut_mesh, &nd_fec);
  auto Project = [](FiniteElementSpace &fespace)
  {
    auto E = std::make_unique<GridFunction>(fespace, true);
    mfem::VectorFunctionCoefficient fr(3, SmoothE), fi(3, JumpE);
    E->Real().ProjectCoefficient(fr);
    E->Imag().ProjectCoefficient(fi);
    E->Real().ExchangeFaceNbrData();  // For the legacy coefficient evaluation
    E->Imag().ExchangeFaceNbrData();
    return E;
  };
  auto E = Project(nd_fespace);
  auto cut_E = Project(cut_nd_fespace);

  const int bdr_attr_max = mesh.Get().bdr_attributes.Max();
  mfem::Array<int> marker(bdr_attr_max), none(bdr_attr_max);
  marker = 0;
  none = 0;
  marker[crack_attr - 1] = 1;
  const double t_i = 2.0e-3, epsilon_i = 10.0;

  // Legacy coefficient integral over the marked boundary elements.
  auto LegacyIntegral = [&](mfem::Coefficient &f, const mfem::ParMesh &pmesh)
  {
    double sum = 0.0;
    for (int be = 0; be < pmesh.GetNBE(); be++)
    {
      if (pmesh.GetBdrAttribute(be) != crack_attr)
      {
        continue;
      }
      auto &T = *const_cast<mfem::ParMesh &>(pmesh).GetBdrElementTransformation(be);
      const auto &ir = mfem::IntRules.Get(T.GetGeometryType(), 2 * order + 2);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
        const auto &ip = ir.IntPoint(q);
        T.SetIntPoint(&ip);
        sum += ip.weight * T.Weight() * f.Eval(T, ip);
      }
    }
    Mpi::GlobalSum(1, &sum, comm);
    return sum;
  };

  for (auto epr_type :
       {InterfaceDielectric::DEFAULT, InterfaceDielectric::MA, InterfaceDielectric::MS})
  {
    CAPTURE(static_cast<int>(epr_type), dielectric);
    SurfaceFunctional sum_sides(mesh, marker, nd_fespace.Get(), mat_op, epr_type, t_i,
                                epsilon_i, &marker);
    SurfaceFunctional average(mesh, marker, nd_fespace.Get(), mat_op, epr_type, t_i,
                              epsilon_i, &none);
    SurfaceFunctional cut_ref(cut_mesh, marker, cut_nd_fespace.Get(), cut_mat_op, epr_type,
                              t_i, epsilon_i);
    REQUIRE(sum_sides.IsValid());
    REQUIRE(cut_ref.IsValid());
    const double val = sum_sides.Eval(*E), ref = cut_ref.Eval(*cut_E),
                 avg = average.Eval(*E);
    const bool qualifies = (epr_type == InterfaceDielectric::DEFAULT) ||
                           ((epr_type == InterfaceDielectric::MA) != dielectric);
    CAPTURE(val, ref, avg);
    if (qualifies)
    {
      CHECK_THAT(val, WithinRel(ref, 1.0e-10));
      CHECK(std::abs(val - avg) > 1.0e-3 * std::abs(val));
    }
    else
    {
      CHECK(std::abs(val) < 1.0e-14);
      CHECK(std::abs(ref) < 1.0e-14);
    }

    // The legacy coefficient path agrees.
    auto Legacy = [&]() -> std::unique_ptr<mfem::Coefficient>
    {
      switch (epr_type)
      {
        case InterfaceDielectric::DEFAULT:
          return std::make_unique<
              InterfaceDielectricCoefficient<InterfaceDielectric::DEFAULT>>(
              *E, mat_op, t_i, epsilon_i, &marker);
        case InterfaceDielectric::MA:
          return std::make_unique<InterfaceDielectricCoefficient<InterfaceDielectric::MA>>(
              *E, mat_op, t_i, epsilon_i, &marker);
        default:
          return std::make_unique<InterfaceDielectricCoefficient<InterfaceDielectric::MS>>(
              *E, mat_op, t_i, epsilon_i, &marker);
      }
    }();
    const double legacy = LegacyIntegral(*Legacy, mesh.Get());
    CHECK_THAT(legacy, WithinRel(val, 1.0e-8) || WithinAbs(val, 1.0e-14));
  }
}

}  // namespace palace
