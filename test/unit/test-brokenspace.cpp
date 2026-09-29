// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <numeric>
#include <random>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "fem/brokenspace.hpp"
#include "fem/coefficient.hpp"
#include "fem/errorindicator.hpp"
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/output_functionals.hpp"
#include "linalg/errorestimator.hpp"
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
// plane x = x_crack, for y <= y_max (the interior boundary has a free edge for y_max < 1).
// For type mfem::Element::PYRAMID, the mesh is that of MakeHexPyramidMesh.
mfem::Mesh MakeCrackedCubeMesh(int n, mfem::Element::Type type, double y_max = 1.0)
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
    bool on_plane = true;
    for (auto v : fv)
    {
      on_plane = on_plane && (std::abs(base.GetVertex(v)[0] - x_crack) < 1.0e-12) &&
                 (base.GetVertex(v)[1] < y_max + 1.0e-12);
    }
    if (on_plane)
    {
      crack_faces.push_back(f);
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

  // Element-wise squared estimates for the gradient and curl flux estimators.
  std::pair<Vector, Vector> Estimate(void (*E)(const mfem::Vector &, mfem::Vector &),
                                     void (*B)(const mfem::Vector &, mfem::Vector &),
                                     const std::vector<int> &crack_attr_list)
  {
    constexpr double tol = 1.0e-14;
    constexpr int max_it = 10000, print = 0;
    constexpr bool use_mg = false;
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

}  // namespace

TEST_CASE("Interior boundary sides", "[brokenspace][Serial][Parallel]")
{
  const auto comm = MPI_COMM_WORLD;
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON,
                             mfem::Element::WEDGE, mfem::Element::PYRAMID,
                             mfem::Element::TRIANGLE, mfem::Element::QUADRILATERAL);
  auto smesh = MakeCrackedCubeMesh(4, type);
  auto part = Partition(smesh, Mpi::Size(comm));
  auto mesh = MakeParMesh(comm, smesh, part);
  const auto &pmesh = mesh.Get();

  const std::vector<int> attr_list = {crack_attr};
  const auto &sides = mesh.GetCrackSides(attr_list);
  REQUIRE(sides.copy.size() == static_cast<std::size_t>(pmesh.GetNE()));
  REQUIRE(sides.split.size() == static_cast<std::size_t>(pmesh.GetNE()));

  // For a plane spanning the domain, exactly one side of every interior boundary entity
  // reads the copy: all elements on one side of the plane touching an entity of the plane
  // have the corresponding bit set and those on the other side do not. The base side of an
  // entity is the one with the smallest global element number, which depends on the
  // partitioning, but the global number of elements with a set bit for a vertex of the
  // plane plus those without must equal the number of elements touching it.
  mfem::Array<int> ev, ee, eo, ef, efo;
  mfem::Vector c(3);
  int num_bad = 0, num_set = 0;
  for (int e = 0; e < pmesh.GetNE(); e++)
  {
    const auto geom = pmesh.GetElementGeometry(e);
    const int nv_e = mfem::Geometry::NumVerts[geom], ne_e = mfem::Geometry::NumEdges[geom];
    pmesh.GetElementVertices(e, ev);
    pmesh.GetElementEdges(e, ee, eo);
    if (pmesh.Dimension() == 3)
    {
      pmesh.GetElementFaces(e, ef, efo);
    }
    else
    {
      ef.SetSize(0);
    }
    std::uint32_t on_plane = 0;
    for (int i = 0; i < ev.Size(); i++)
    {
      if (std::abs(pmesh.GetVertex(ev[i])[0] - x_crack) < 1.0e-12)
      {
        on_plane |= (1u << i);
      }
    }
    mfem::Array<int> verts;
    for (int i = 0; i < ee.Size(); i++)
    {
      pmesh.GetEdgeVertices(ee[i], verts);
      if (std::abs(pmesh.GetVertex(verts[0])[0] - x_crack) < 1.0e-12 &&
          std::abs(pmesh.GetVertex(verts[1])[0] - x_crack) < 1.0e-12)
      {
        on_plane |= (1u << (nv_e + i));
      }
    }
    for (int i = 0; i < ef.Size(); i++)
    {
      pmesh.GetFaceVertices(ef[i], verts);
      bool face_on_plane = true;
      for (auto v : verts)
      {
        face_on_plane =
            face_on_plane && (std::abs(pmesh.GetVertex(v)[0] - x_crack) < 1.0e-12);
      }
      if (face_on_plane)
      {
        on_plane |= (1u << (nv_e + ne_e + i));
      }
    }
    // Bits are only set for entities on the plane (all of which are split).
    num_bad += ((sides.copy[e] & ~on_plane) != 0) + (sides.split[e] != on_plane) +
               ((sides.copy[e] & ~sides.split[e]) != 0);
    num_set += (sides.copy[e] != 0);
  }
  Mpi::GlobalSum(1, &num_bad, comm);
  Mpi::GlobalSum(1, &num_set, comm);
  CHECK(num_bad == 0);
  CHECK(num_set > 0);

  // No interior boundary: no sides.
  const std::vector<int> other_list = {1};
  const auto &no_sides = mesh.GetCrackSides(other_list);
  int num_nonzero = 0;
  for (std::size_t e = 0; e < no_sides.copy.size(); e++)
  {
    num_nonzero += (no_sides.copy[e] != 0) + (no_sides.split[e] != 0);
  }
  Mpi::GlobalSum(1, &num_nonzero, comm);
  CHECK(num_nonzero == 0);
}

TEST_CASE("Hanging entities", "[brokenspace][Serial][Parallel]")
{
  const auto comm = MPI_COMM_WORLD;
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON,
                             mfem::Element::TRIANGLE, mfem::Element::QUADRILATERAL);
  auto smesh = MakeCrackedCubeMesh(4, type);
  smesh.EnsureNCMesh(true);
  auto part = Partition(smesh, Mpi::Size(comm));
  mfem::ParMesh pmesh(comm, smesh, part.data());
  CHECK(!mesh::HasHangingEntities(pmesh));

  // Uniform refinement leaves no hanging entities, local refinement does.
  mfem::Array<int> all(pmesh.GetNE());
  std::iota(all.begin(), all.end(), 0);
  pmesh.GeneralRefinement(all, 1, 0);
  CHECK(!mesh::HasHangingEntities(pmesh));
  mfem::Array<int> first;
  if (Mpi::Root(comm))
  {
    first.Append(0);
  }
  pmesh.GeneralRefinement(first, 1, 0);
  CHECK(mesh::HasHangingEntities(pmesh));
}

TEST_CASE("Interior boundary sides through refinement", "[brokenspace][Serial][Parallel]")
{
  // After a uniform nonconforming refinement, which leaves no hanging entities, the sides
  // inherited through the refinement are those computed on the refined mesh (up to the
  // choice of the base side). The mesh has no nodes, since MFEM does not support the
  // refinement of the nodes of a mesh with pyramids.
  const auto comm = MPI_COMM_WORLD;
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON,
                             mfem::Element::WEDGE, mfem::Element::PYRAMID,
                             mfem::Element::TRIANGLE, mfem::Element::QUADRILATERAL);
  const double y_max = GENERATE(1.0, 0.5);
  auto smesh = MakeCrackedCubeMesh(4, type, y_max);
  smesh.EnsureNCMesh(true);
  auto part = Partition(smesh, Mpi::Size(comm));
  mfem::ParMesh pmesh(comm, smesh, part.data());
  mfem::Array<int> marker(crack_attr);
  marker = 0;
  marker[crack_attr - 1] = 1;
  const auto coarse_sides = mesh::ComputeCrackSides(pmesh, marker);
  mfem::Array<int> all(pmesh.GetNE());
  std::iota(all.begin(), all.end(), 0);
  pmesh.GeneralRefinement(all, 1, 0);
  REQUIRE(!mesh::HasHangingEntities(pmesh));
  const auto inherited =
      mesh::InheritCrackSides(pmesh, pmesh.GetRefinementTransforms(), coarse_sides);
  const auto computed = mesh::ComputeCrackSides(pmesh, marker);

  // The split entities agree, and so do the copies of split entities, on every element or
  // on none (the base side is the other one). The inherited sides may also read copies of
  // entities which are not split in the closure of split entities of the parent, which has
  // no effect when they are not constrained (there are no hanging entities).
  int num_split = 0, num_diff_split = 0, num_same_copy = 0, num_flipped_copy = 0;
  for (int e = 0; e < pmesh.GetNE(); e++)
  {
    const auto split = computed.split[e];
    num_split += (split != 0);
    num_diff_split += (inherited.split[e] != split);
    if (split)
    {
      const auto copy = inherited.copy[e] & split;
      num_same_copy += (copy == computed.copy[e]);
      num_flipped_copy += (copy == (split & ~computed.copy[e]));
    }
  }
  Mpi::GlobalSum(1, &num_split, comm);
  Mpi::GlobalSum(1, &num_diff_split, comm);
  Mpi::GlobalSum(1, &num_same_copy, comm);
  Mpi::GlobalSum(1, &num_flipped_copy, comm);
  INFO("type " << type << ", y_max " << y_max << ": split " << num_split << ", diff. split "
               << num_diff_split << ", same copy " << num_same_copy << ", flipped copy "
               << num_flipped_copy);
  CHECK(num_split > 0);
  CHECK(num_diff_split == 0);
  CHECK((num_same_copy == num_split || num_flipped_copy == num_split));
}

TEST_CASE("Broken space prolongation", "[brokenspace][Serial][Parallel]")
{
  const auto comm = MPI_COMM_WORLD;
  const int order = GENERATE(1, 2);
  auto smesh = MakeCrackedCubeMesh(4, mfem::Element::TETRAHEDRON);
  auto part = Partition(smesh, Mpi::Size(comm));
  auto mesh = MakeParMesh(comm, smesh, part);
  mfem::ND_FECollection nd_fec(order, 3);
  mfem::RT_FECollection rt_fec(order - 1, 3);
  mfem::H1_FECollection h1_fec(order, 3);
  const std::vector<int> attr_list = {crack_attr};
  for (const mfem::FiniteElementCollection *fec :
       std::vector<const mfem::FiniteElementCollection *>{&nd_fec, &rt_fec, &h1_fec})
  {
    FiniteElementSpace base_fespace(mesh, fec);
    const auto tsize0 = base_fespace.GlobalTrueVSize();
    FiniteElementSpace fespace(base_fespace, mesh.GetCrackSides(attr_list));
    REQUIRE(fespace.IsBroken());
    CHECK(fespace.GlobalTrueVSize() > tsize0);

    // The broken view shares the MFEM space, and leaves the given space unchanged.
    CHECK(&fespace.Get() == &base_fespace.Get());
    CHECK(!base_fespace.IsBroken());
    CHECK(base_fespace.GlobalTrueVSize() == tsize0);
    CHECK(base_fespace.GetProlongationMatrix() != fespace.GetProlongationMatrix());

    // Compare with the global size of the space on the cut mesh.
    auto cut = CutMesh(smesh);
    auto cut_mesh = MakeParMesh(comm, cut, part);
    FiniteElementSpace cut_fespace(cut_mesh, fec);
    CHECK(fespace.GlobalTrueVSize() == cut_fespace.GlobalTrueVSize());

    // Adjointness of the prolongation and its transpose.
    const Operator &P = *fespace.GetProlongationMatrix();
    Vector x(P.Width()), y(P.Height()), Px(P.Height()), Pty(P.Width());
    std::mt19937 gen(Mpi::Rank(comm) + 1);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (int i = 0; i < x.Size(); i++)
    {
      x(i) = dist(gen);
    }
    for (int i = 0; i < y.Size(); i++)
    {
      y(i) = dist(gen);
    }
    P.Mult(x, Px);
    P.MultTranspose(y, Pty);
    // The L-vector inner product is not a global inner product (shared L-DOFs are counted
    // on each process), but y · P x = Pᵀ y · x holds summed over processes.
    double dots[2] = {y * Px, Pty * x};
    Mpi::GlobalSum(2, dots, comm);
    CHECK_THAT(dots[0], WithinRel(dots[1], 1.0e-12));
  }
}

TEST_CASE("Broken space error estimators",
          "[brokenspace][errorestimator][Serial][Parallel]")
{
  const auto comm = MPI_COMM_WORLD;
  const int order = GENERATE(1, 2);
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON,
                             mfem::Element::PYRAMID);
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
  mesh.RefineCrackSides();
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
  const int order = GENERATE(1, 2);
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
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
  uncut_setup.mesh.GetCrackSides(attr_list);

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
    REQUIRE(mesh::HasHangingEntities(uncut_setup.mesh.Get()));
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
  }

  SECTION("Exact recovery of a jump with hanging entities across the interior boundary")
  {
    // Refinement of one side only, then of part of the other side.
    RefineWhere(uncut_setup.mesh, [](const mfem::Vector &c)
                { return c(0) < x_crack && c(0) > 0.25 && c(2) < 0.75; });
    UpdateSetup(uncut_setup.mesh, uncut_setup.nd_fespaces, uncut_setup.rt_fespaces);
    uncut_setup.mesh.GetCrackSides(attr_list);
    RefineWhere(uncut_setup.mesh, [](const mfem::Vector &c)
                { return c(0) > x_crack && c(0) < 0.75 && c(1) < 0.4; });
    UpdateSetup(uncut_setup.mesh, uncut_setup.nd_fespaces, uncut_setup.rt_fespaces);
    REQUIRE(mesh::HasHangingEntities(uncut_setup.mesh.Get()));
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
      other_setup.mesh.GetCrackSides(attr_list);
      for (auto *setup : {&uncut_setup, &other_setup})
      {
        RefineWhere(setup->mesh, [](const mfem::Vector &c)
                    { return c(0) < x_crack && c(0) > 0.25 && c(2) < 0.75; });
        UpdateSetup(setup->mesh, setup->nd_fespaces, setup->rt_fespaces);
        setup->mesh.GetCrackSides(attr_list);
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

TEST_CASE("Broken space error estimators (partial interior boundary)",
          "[brokenspace][errorestimator][Serial][Parallel]")
{
  // An interior boundary with a free edge (at y = 0.5), where the recovery remains
  // continuous, and a partitioning where some processes do not touch the interior boundary.
  const auto comm = MPI_COMM_WORLD;
  const int order = GENERATE(1, 2);
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
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
    uncut_setup.mesh.GetCrackSides(attr_list);
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
  const int order = GENERATE(1, 2);
  const auto type = GENERATE(mfem::Element::TRIANGLE, mfem::Element::QUADRILATERAL);
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
