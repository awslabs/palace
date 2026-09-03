// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <mfem.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include "fem/boundary_derived_field_bundle.hpp"
#include "fem/boundary_physical_trace.hpp"
#include "fem/coefficient.hpp"
#include "fem/face_sampling_plan.hpp"
#include "fem/facenbrexchange.hpp"
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/output_functionals.hpp"
#include "fem/point_field_evaluator.hpp"
#include "fixtures.hpp"
#include "models/boundarymodeoperator.hpp"
#include "models/domainpostoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/postoperator.hpp"
#include "models/strattonchu.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/geodata.hpp"
#include "utils/labels.hpp"
#include "utils/units.hpp"

namespace palace
{

namespace
{

std::unique_ptr<Mesh> LoadTestMesh(MPI_Comm comm, const std::string &file)
{
  mfem::Mesh smesh(std::string(PALACE_TEST_DATA_DIR) + "/mesh/" + file, 1, 1);
  smesh.EnsureNodes();
  REQUIRE(Mpi::Size(comm) <= smesh.GetNE());
  auto pmesh = std::make_unique<mfem::ParMesh>(comm, smesh);
  return std::make_unique<Mesh>(std::move(pmesh));
}

// Reference (host, mfem-based) computation of the surface area with the same
// quadrature order as SurfaceFunctional (mesh::GetSurfaceArea uses a different
// quadrature order, which only matters for curved boundary elements).
double RefSurfaceArea(mfem::ParMesh &pmesh, const mfem::Array<int> &marker)
{
  double area = 0.0;
  for (int i = 0; i < pmesh.GetNBE(); i++)
  {
    if (!marker[pmesh.GetBdrAttribute(i) - 1])
    {
      continue;
    }
    auto *T = pmesh.GetBdrElementTransformation(i);
    const auto &ir = mfem::IntRules.Get(pmesh.GetBdrElementGeometry(i),
                                        fem::DefaultIntegrationOrder::Get(*T));
    for (int q = 0; q < ir.GetNPoints(); q++)
    {
      const auto &ip = ir.IntPoint(q);
      T->SetIntPoint(&ip);
      area += ip.weight * T->Weight();
    }
  }
  Mpi::GlobalSum(1, &area, pmesh.GetComm());
  return area;
}

// Reference (host, mfem-based) computation of ∫ |E|² dS over the marked boundary
// elements, evaluating E from the attached volume element at each quadrature point in
// the same way as the legacy BdrGridFunctionCoefficient-based postprocessing.
double RefSurfaceHCurlNorm2(mfem::ParMesh &pmesh, const mfem::ParGridFunction &E,
                            const mfem::Array<int> &marker)
{
  double sum = 0.0;
  mfem::FaceElementTransformations FET;
  mfem::IsoparametricTransformation T1, T2;
  mfem::Vector Ev(pmesh.SpaceDimension());
  for (int i = 0; i < pmesh.GetNBE(); i++)
  {
    if (!marker[pmesh.GetBdrAttribute(i) - 1])
    {
      continue;
    }
    auto *T = pmesh.GetBdrElementTransformation(i);
    const auto &ir = mfem::IntRules.Get(pmesh.GetBdrElementGeometry(i),
                                        fem::DefaultIntegrationOrder::Get(*T));
    for (int q = 0; q < ir.GetNPoints(); q++)
    {
      const auto &ip = ir.IntPoint(q);
      T->SetIntPoint(&ip);
      BdrGridFunctionCoefficient::GetBdrElementNeighborTransformations(i, pmesh, FET, T1,
                                                                       T2, &ip);
      E.GetVectorValue(*FET.Elem1, FET.Elem1->GetIntPoint(), Ev);
      sum += ip.weight * T->Weight() * (Ev * Ev);
    }
  }
  Mpi::GlobalSum(1, &sum, pmesh.GetComm());
  return sum;
}

// Reference (host, mfem-based) integration of a legacy boundary coefficient over the
// marked boundary elements, mirroring SurfacePostOperator::GetLocalSurfaceIntegral.
double RefSurfaceCoefficientIntegral(mfem::ParMesh &pmesh, mfem::Coefficient &f,
                                     const mfem::Array<int> &marker)
{
  double sum = 0.0;
  for (int i = 0; i < pmesh.GetNBE(); i++)
  {
    if (!marker[pmesh.GetBdrAttribute(i) - 1])
    {
      continue;
    }
    auto *T = pmesh.GetBdrElementTransformation(i);
    const auto &ir = mfem::IntRules.Get(pmesh.GetBdrElementGeometry(i),
                                        fem::DefaultIntegrationOrder::Get(*T));
    for (int q = 0; q < ir.GetNPoints(); q++)
    {
      const auto &ip = ir.IntPoint(q);
      T->SetIntPoint(&ip);
      sum += ip.weight * T->Weight() * f.Eval(*T, ip);
    }
  }
  Mpi::GlobalSum(1, &sum, pmesh.GetComm());
  return sum;
}

// Build a 2D mesh with two material regions (attribute 1: y < 0.5, vacuum; attribute
// 2: y >= 0.5, dielectric) and interior boundary elements (attribute 7) added on the
// material interface at y = 0.5.
std::unique_ptr<Mesh> MakeInterfaceMesh2D(MPI_Comm comm, mfem::Element::Type elem_type)
{
  mfem::Mesh smesh = mfem::Mesh::MakeCartesian2D(2, 2, elem_type, false, 1.0, 1.0);
  for (int e = 0; e < smesh.GetNE(); e++)
  {
    mfem::Vector center(2);
    smesh.GetElementCenter(e, center);
    smesh.SetAttribute(e, (center(1) < 0.5) ? 1 : 2);
  }
  // Add interior boundary elements on the material interface.
  for (int f = 0; f < smesh.GetNumFaces(); f++)
  {
    int e1, e2;
    smesh.GetFaceElements(f, &e1, &e2);
    if (e1 >= 0 && e2 >= 0 && smesh.GetAttribute(e1) != smesh.GetAttribute(e2))
    {
      auto *face_elem = smesh.GetFace(f)->Duplicate(&smesh);
      face_elem->SetAttribute(7);
      smesh.AddBdrElement(face_elem);
    }
  }
  smesh.FinalizeTopology();
  smesh.Finalize();
  smesh.SetAttributes();
  smesh.EnsureNodes();
  REQUIRE(Mpi::Size(comm) <= smesh.GetNE());
  auto pmesh = std::make_unique<mfem::ParMesh>(comm, smesh);
  return std::make_unique<Mesh>(std::move(pmesh));
}

// Two-dimensional counterpart of MakeNCInterfaceMesh. Refining only the lower material
// region leaves interface subfaces attached to a coarse upper-material neighbor,
// including depth-2 reference maps when requested.
std::unique_ptr<Mesh> MakeNCInterfaceMesh2D(MPI_Comm comm, mfem::Element::Type elem_type,
                                            int refinement_depth, int coarse_size = 2)
{
  mfem::Mesh smesh =
      mfem::Mesh::MakeCartesian2D(coarse_size, coarse_size, elem_type, false, 1.0, 1.0);
  for (int e = 0; e < smesh.GetNE(); e++)
  {
    mfem::Vector center(2);
    smesh.GetElementCenter(e, center);
    smesh.SetAttribute(e, (center(1) < 0.5) ? 1 : 2);
  }
  for (int f = 0; f < smesh.GetNumFaces(); f++)
  {
    int e1, e2;
    smesh.GetFaceElements(f, &e1, &e2);
    if (e1 >= 0 && e2 >= 0 && smesh.GetAttribute(e1) != smesh.GetAttribute(e2))
    {
      auto *face_elem = smesh.GetFace(f)->Duplicate(&smesh);
      face_elem->SetAttribute(7);
      smesh.AddBdrElement(face_elem);
    }
  }
  smesh.FinalizeTopology();
  smesh.Finalize();
  smesh.SetAttributes();
  smesh.EnsureNodes();
  smesh.EnsureNCMesh(true);

  REQUIRE(refinement_depth >= 1);
  for (int level = 0; level < refinement_depth; level++)
  {
    mfem::Array<int> refs;
    for (int e = 0; e < smesh.GetNE(); e++)
    {
      mfem::Vector center(2);
      smesh.GetElementCenter(e, center);
      if (center(1) < 0.5)
      {
        refs.Append(e);
      }
    }
    smesh.GeneralRefinement(refs, 1);
  }

  REQUIRE(Mpi::Size(comm) <= smesh.GetNE());
  auto pmesh = std::make_unique<mfem::ParMesh>(comm, smesh);
  return std::make_unique<Mesh>(std::move(pmesh));
}

class TestBoundaryModePostOperator : public PostOperator<ProblemType::BOUNDARYMODE>
{
public:
  using PostOperator<ProblemType::BOUNDARYMODE>::PostOperator;
  const Measurement &Cache() const { return measurement_cache; }
};

// Build a 3D mesh with two material regions (attribute 1: z < 0.5, vacuum; attribute
// 2: z >= 0.5, dielectric) and interior boundary elements (attribute 7) added on the
// material interface at z = 0.5.
std::unique_ptr<Mesh> MakeInterfaceMesh(MPI_Comm comm, mfem::Element::Type elem_type)
{
  mfem::Mesh smesh = mfem::Mesh::MakeCartesian3D(2, 2, 2, elem_type);
  for (int e = 0; e < smesh.GetNE(); e++)
  {
    mfem::Vector center(3);
    smesh.GetElementCenter(e, center);
    smesh.SetAttribute(e, (center(2) < 0.5) ? 1 : 2);
  }
  // Add interior boundary elements on the material interface.
  for (int f = 0; f < smesh.GetNumFaces(); f++)
  {
    int e1, e2;
    smesh.GetFaceElements(f, &e1, &e2);
    if (e1 >= 0 && e2 >= 0 && smesh.GetAttribute(e1) != smesh.GetAttribute(e2))
    {
      auto *face_elem = smesh.GetFace(f)->Duplicate(&smesh);
      face_elem->SetAttribute(7);
      smesh.AddBdrElement(face_elem);
    }
  }
  smesh.FinalizeTopology();
  smesh.Finalize();
  smesh.SetAttributes();
  smesh.EnsureNodes();
  REQUIRE(Mpi::Size(comm) <= smesh.GetNE());
  auto pmesh = std::make_unique<mfem::ParMesh>(comm, smesh);
  return std::make_unique<Mesh>(std::move(pmesh));
}

// Build a 3D mesh with two material regions and interior boundary elements (attribute
// 7) on the material interface at z = 0.5, then apply nonconformal refinement to a
// subset of elements so that marked surfaces (exterior and interior) contain
// nonconformal (master/slave) faces.
std::unique_ptr<Mesh> MakeNCInterfaceMesh(MPI_Comm comm, mfem::Element::Type elem_type,
                                          int refinement_depth = 1)
{
  mfem::Mesh smesh = mfem::Mesh::MakeCartesian3D(2, 2, 2, elem_type);
  for (int e = 0; e < smesh.GetNE(); e++)
  {
    mfem::Vector center(3);
    smesh.GetElementCenter(e, center);
    smesh.SetAttribute(e, (center(2) < 0.5) ? 1 : 2);
  }
  for (int f = 0; f < smesh.GetNumFaces(); f++)
  {
    int e1, e2;
    smesh.GetFaceElements(f, &e1, &e2);
    if (e1 >= 0 && e2 >= 0 && smesh.GetAttribute(e1) != smesh.GetAttribute(e2))
    {
      auto *face_elem = smesh.GetFace(f)->Duplicate(&smesh);
      face_elem->SetAttribute(7);
      smesh.AddBdrElement(face_elem);
    }
  }
  smesh.FinalizeTopology();
  smesh.Finalize();
  smesh.SetAttributes();
  smesh.EnsureNodes();
  smesh.EnsureNCMesh(true);

  // Nonconformally refine the elements in the x < 0.5 half (crosses both the interior
  // interface and the exterior boundaries), creating master/slave faces on the marked
  // surfaces. Repeating this fixed selection produces a deterministic depth-2 case
  // without changing the default depth-1 fixtures used by the broader test matrix.
  REQUIRE(refinement_depth >= 1);
  for (int level = 0; level < refinement_depth; level++)
  {
    mfem::Array<int> refs;
    for (int e = 0; e < smesh.GetNE(); e++)
    {
      mfem::Vector center(3);
      smesh.GetElementCenter(e, center);
      if (center(0) < 0.5)
      {
        refs.Append(e);
      }
    }
    smesh.GeneralRefinement(refs, 1);
  }

  REQUIRE(Mpi::Size(comm) <= smesh.GetNE());
  auto pmesh = std::make_unique<mfem::ParMesh>(comm, smesh);
  return std::make_unique<Mesh>(std::move(pmesh));
}

// Refine only one material side so the retained interior boundary is represented by
// multiple fine leaves attached to one coarse volume owner. This specifically exercises
// the 3D coarse-face union rather than merely putting NC volume faces near a surface.
std::unique_ptr<Mesh> MakeNCCoarseInterfaceMesh3D(MPI_Comm comm,
                                                  mfem::Element::Type elem_type,
                                                  int refinement_depth, int coarse_size = 2)
{
  mfem::Mesh smesh = mfem::Mesh::MakeCartesian3D(coarse_size, coarse_size, 2, elem_type);
  for (int e = 0; e < smesh.GetNE(); e++)
  {
    mfem::Vector center(3);
    smesh.GetElementCenter(e, center);
    smesh.SetAttribute(e, (center(2) < 0.5) ? 1 : 2);
  }
  for (int f = 0; f < smesh.GetNumFaces(); f++)
  {
    int e1, e2;
    smesh.GetFaceElements(f, &e1, &e2);
    if (e1 >= 0 && e2 >= 0 && smesh.GetAttribute(e1) != smesh.GetAttribute(e2))
    {
      auto *face_elem = smesh.GetFace(f)->Duplicate(&smesh);
      face_elem->SetAttribute(7);
      smesh.AddBdrElement(face_elem);
    }
  }
  smesh.FinalizeTopology();
  smesh.Finalize();
  smesh.SetAttributes();
  smesh.EnsureNodes();
  smesh.EnsureNCMesh(true);

  REQUIRE(refinement_depth >= 1);
  for (int level = 0; level < refinement_depth; level++)
  {
    mfem::Array<int> refs;
    for (int e = 0; e < smesh.GetNE(); e++)
    {
      mfem::Vector center(3);
      smesh.GetElementCenter(e, center);
      if (center(2) < 0.5)
      {
        refs.Append(e);
      }
    }
    smesh.GeneralRefinement(refs, 1);
  }

  REQUIRE(Mpi::Size(comm) <= smesh.GetNE());
  auto pmesh = std::make_unique<mfem::ParMesh>(comm, smesh);
  return std::make_unique<Mesh>(std::move(pmesh));
}

enum class MixedCoarseInterface
{
  PRISM_QUAD,
  PRISM_TRIANGLE,
  PYRAMID_QUAD,
  PYRAMID_TRIANGLE
};

mfem::Geometry::Type CoarseGeometry(MixedCoarseInterface interface)
{
  return interface == MixedCoarseInterface::PRISM_QUAD ||
                 interface == MixedCoarseInterface::PRISM_TRIANGLE
             ? mfem::Geometry::PRISM
             : mfem::Geometry::PYRAMID;
}

mfem::Geometry::Type CoarseFaceGeometry(MixedCoarseInterface interface)
{
  return interface == MixedCoarseInterface::PRISM_QUAD ||
                 interface == MixedCoarseInterface::PYRAMID_QUAD
             ? mfem::Geometry::SQUARE
             : mfem::Geometry::TRIANGLE;
}

// Add retained interior boundaries to replicated connected mixed-element zoos, then
// refine one neighbor so a wedge or pyramid owns coarse NC leaves on both supported
// triangular and quadrilateral facet families.
std::unique_ptr<Mesh> MakeNCMixedCoarseOwnerMesh(MPI_Comm comm,
                                                 MixedCoarseInterface interface,
                                                 int refinement_depth)
{
  mfem::Mesh source(std::string(PALACE_TEST_DATA_DIR) + "/mesh/tinyzoo-3d.mesh", 1, 1);
  REQUIRE(source.GetNE() == 4);
  // Replicate disconnected zoo components so MPI partitions with more than two ranks do
  // not leave a rank with an under-populated mixed NC component during Mesh setup.
  const int copies = std::max(1, Mpi::Size(comm));
  const int vertices_per_copy = source.GetNV() + 1;
  const int elements_per_copy = source.GetNE() + 1;
  mfem::Mesh smesh(3, copies * vertices_per_copy, copies * elements_per_copy, 0, 3);
  mfem::Array<int> vertices;
  for (int copy = 0; copy < copies; copy++)
  {
    const int vertex_offset = copy * vertices_per_copy;
    for (int v = 0; v < source.GetNV(); v++)
    {
      const double *x = source.GetVertex(v);
      smesh.AddVertex(x[0] + 3.0 * copy, x[1], x[2]);
    }
    // Attach one extra tetrahedron below the wedge's z=0 triangular face {4,5,1}.
    // This makes that wedge facet an interior interface whose tetrahedral neighbor can
    // be refined while the wedge remains coarse.
    const int extra_vertex = vertex_offset + source.GetNV();
    smesh.AddVertex(1.25 + 3.0 * copy, 0.5, -1.0);
    for (int e = 0; e < source.GetNE(); e++)
    {
      source.GetElementVertices(e, vertices);
      for (int j = 0; j < vertices.Size(); j++)
      {
        vertices[j] += vertex_offset;
      }
      const int attr = e + 1;
      switch (source.GetElementBaseGeometry(e))
      {
        case mfem::Geometry::TETRAHEDRON:
          smesh.AddTet(vertices.GetData(), attr);
          break;
        case mfem::Geometry::CUBE:
          smesh.AddHex(vertices.GetData(), attr);
          break;
        case mfem::Geometry::PRISM:
          smesh.AddWedge(vertices.GetData(), attr);
          break;
        case mfem::Geometry::PYRAMID:
          smesh.AddPyramid(vertices.GetData(), attr);
          break;
        default:
          MFEM_ABORT("Unexpected geometry in mixed NC zoo fixture!");
      }
    }
    const int prism_triangle_neighbor[4] = {vertex_offset + 4, vertex_offset + 5,
                                            vertex_offset + 1, extra_vertex};
    smesh.AddTet(prism_triangle_neighbor, 5);
  }
  smesh.FinalizeTopology();
  smesh.Finalize();
  // Retain every mixed interior face as a point-output boundary before NC refinement.
  for (int f = 0; f < smesh.GetNumFaces(); f++)
  {
    int e1, e2;
    smesh.GetFaceElements(f, &e1, &e2);
    if (e1 >= 0 && e2 >= 0)
    {
      auto *face_elem = smesh.GetFace(f)->Duplicate(&smesh);
      face_elem->SetAttribute(7);
      smesh.AddBdrElement(face_elem);
    }
  }
  smesh.FinalizeTopology();
  smesh.Finalize();
  smesh.SetAttributes();
  smesh.EnsureNodes();
  smesh.EnsureNCMesh(true);

  // Attributes are hex=1, wedge=2, pyramid=3, original tet=4, added tet=5.
  const int refined_attribute = interface == MixedCoarseInterface::PRISM_QUAD       ? 1
                                : interface == MixedCoarseInterface::PRISM_TRIANGLE ? 5
                                : interface == MixedCoarseInterface::PYRAMID_QUAD   ? 2
                                                                                    : 4;
  REQUIRE(refinement_depth >= 1);
  for (int level = 0; level < refinement_depth; level++)
  {
    mfem::Array<int> refs;
    for (int e = 0; e < smesh.GetNE(); e++)
    {
      if (smesh.GetAttribute(e) == refined_attribute)
      {
        refs.Append(e);
      }
    }
    REQUIRE(refs.Size() > 0);
    smesh.GeneralRefinement(refs, 1);
  }

  REQUIRE(Mpi::Size(comm) <= smesh.GetNE());
  std::vector<int> partitioning;
  if (Mpi::Size(comm) > 1)
  {
    // Force every targeted mixed coarse/fine interface across ranks 0 and 1. Distribute
    // unrelated elements over the remaining ranks so MPI-8 has no empty partition.
    const int coarse_attribute = CoarseGeometry(interface) == mfem::Geometry::PRISM ? 2 : 3;
    partitioning.resize(smesh.GetNE());
    for (int e = 0; e < smesh.GetNE(); e++)
    {
      const int attr = smesh.GetAttribute(e);
      if (attr == coarse_attribute)
      {
        partitioning[e] = 0;
      }
      else if (attr == refined_attribute)
      {
        partitioning[e] = 1;
      }
      else if (Mpi::Size(comm) == 2)
      {
        partitioning[e] = e % 2;
      }
      else
      {
        partitioning[e] = 2 + e % (Mpi::Size(comm) - 2);
      }
    }
  }
  auto pmesh = std::make_unique<mfem::ParMesh>(
      comm, smesh, partitioning.empty() ? nullptr : partitioning.data());
  return std::make_unique<Mesh>(std::move(pmesh));
}

void CheckBoundaryTraceNCUnion(
    std::unique_ptr<Mesh> mesh,
    mfem::Geometry::Type expected_coarse_geometry = mfem::Geometry::INVALID,
    mfem::Geometry::Type expected_coarse_face_geometry = mfem::Geometry::INVALID,
    bool expect_ghost_union = false)
{
  MPI_Comm comm = mesh->GetComm();
  auto &pmesh = mesh->Get();
  REQUIRE(pmesh.Nonconforming());
  const int dim = pmesh.Dimension();
  mfem::ND_FECollection nd_fec(2, dim);
  FiniteElementSpace nd_fespace(*mesh, &nd_fec);
  GridFunction E(nd_fespace);
  mfem::VectorFunctionCoefficient field(dim,
                                        [](const mfem::Vector &x, mfem::Vector &v)
                                        {
                                          v(0) = std::sin(x(1)) + 0.25 * x(0);
                                          v(1) = std::cos(x(0)) + x(0) * x(1);
                                          if (v.Size() == 3)
                                          {
                                            v(2) = std::sin(x(0) + x(1)) + 0.5 * x(2);
                                          }
                                        });
  E.Real().ProjectCoefficient(field);
  E.Real().ExchangeFaceNbrData();

  mfem::Array<int> marker(pmesh.bdr_attributes.Max());
  marker = 1;
  constexpr int lod = 2;
  auto sampling_plan = std::make_shared<FaceSamplingPlan>(*mesh, marker, lod);
  long long expected_leaf_points = 0, expected_union_points = 0;
  long long expected_routes = 0, expected_requests_before = 0, expected_requests_after = 0;
  long long expected_local_groups = 0, expected_ghost_groups = 0;
  std::set<std::tuple<bool, int, int>> expected_routing_patterns;
  for (const auto &group : sampling_plan->UnionGroups())
  {
    expected_union_points += static_cast<long long>(group.points.size());
    expected_local_groups += !group.ghost;
    expected_ghost_groups += group.ghost;
    for (const auto &route : group.routes)
    {
      expected_leaf_points += static_cast<long long>(route.canonical_to_union.size());
      expected_routes += static_cast<long long>(route.canonical_to_union.size());
      expected_routing_patterns.emplace(group.ghost, group.side,
                                        static_cast<int>(route.canonical_to_union.size()));
      if (group.ghost)
      {
        expected_requests_before++;
      }
    }
    expected_requests_after += group.ghost;
  }
  bool has_union = !sampling_plan->UnionGroups().empty();
  Mpi::GlobalOr(1, &has_union, comm);
  REQUIRE(has_union);
  if (expected_coarse_geometry != mfem::Geometry::INVALID)
  {
    bool has_expected_coarse_union = std::any_of(
        sampling_plan->UnionGroups().begin(), sampling_plan->UnionGroups().end(),
        [&](const auto &group)
        {
          if (group.geometry != expected_coarse_geometry)
          {
            return false;
          }
          return expected_coarse_face_geometry == mfem::Geometry::INVALID ||
                 std::all_of(group.routes.begin(), group.routes.end(),
                             [&](const auto &route)
                             {
                               return sampling_plan->Entries()[route.entry].bdr_geom ==
                                      expected_coarse_face_geometry;
                             });
        });
    Mpi::GlobalOr(1, &has_expected_coarse_union, comm);
    REQUIRE(has_expected_coarse_union);
    if (expect_ghost_union)
    {
      bool has_expected_ghost_union = std::any_of(
          sampling_plan->UnionGroups().begin(), sampling_plan->UnionGroups().end(),
          [&](const auto &group)
          {
            if (!group.ghost || group.geometry != expected_coarse_geometry)
            {
              return false;
            }
            return expected_coarse_face_geometry == mfem::Geometry::INVALID ||
                   std::all_of(group.routes.begin(), group.routes.end(),
                               [&](const auto &route)
                               {
                                 return sampling_plan->Entries()[route.entry].bdr_geom ==
                                        expected_coarse_face_geometry;
                               });
          });
      Mpi::GlobalOr(1, &has_expected_ghost_union, comm);
      REQUIRE(has_expected_ghost_union);
    }
  }
  REQUIRE((expected_leaf_points > expected_union_points || Mpi::Size(comm) > 1));

  const long long leaf_before = BoundaryPhysicalTraceCache::LeafPointCount();
  const long long union_before = BoundaryPhysicalTraceCache::UnionPointCount();
  const long long duplicate_before = BoundaryPhysicalTraceCache::DuplicatePointCount();
  const long long local_before = BoundaryPhysicalTraceCache::LocalUnionGroupCount();
  const long long local_operators_before =
      BoundaryPhysicalTraceCache::LocalUnionOperatorCount();
  const long long ghost_before = BoundaryPhysicalTraceCache::GhostUnionGroupCount();
  const long long requests_before = BoundaryPhysicalTraceCache::RequestCountBefore();
  const long long requests_after = BoundaryPhysicalTraceCache::RequestCountAfter();
  const long long routes_before = BoundaryPhysicalTraceCache::UnionRouteCount();
  const long long routing_operators_before =
      BoundaryPhysicalTraceCache::RoutingOperatorCount();
  auto trace_cache =
      std::make_shared<BoundaryPhysicalTraceCache>(*mesh, marker, sampling_plan);
  PointFieldEvaluator legacy(PointFieldEvaluator::Kind::FIELD_E, *mesh, marker, nd_fespace,
                             lod, sampling_plan);
  PointFieldEvaluator traced(PointFieldEvaluator::Kind::FIELD_E, *mesh, marker, nd_fespace,
                             lod, sampling_plan, trace_cache);
  REQUIRE(legacy.IsValid());
  REQUIRE(traced.IsValid());
  Vector expected(legacy.BufferSize()), actual(traced.BufferSize());
  expected.UseDevice(true);
  actual.UseDevice(true);
  legacy.EvalBuffer(E.Real(), expected);
  traced.EvalBuffer(E.Real(), actual);
  REQUIRE(expected.Size() == actual.Size());
  const double *expected_values = expected.HostRead();
  const double *actual_values = actual.HostRead();
  for (int i = 0; i < expected.Size(); i++)
  {
    CAPTURE(i, expected_values[i], actual_values[i]);
    CHECK(actual_values[i] ==
          Catch::Approx(expected_values[i]).epsilon(1.0e-10).margin(1.0e-13));
  }

  CHECK(BoundaryPhysicalTraceCache::LeafPointCount() == leaf_before + expected_leaf_points);
  CHECK(BoundaryPhysicalTraceCache::UnionPointCount() ==
        union_before + expected_union_points);
  CHECK(BoundaryPhysicalTraceCache::DuplicatePointCount() ==
        duplicate_before + expected_leaf_points - expected_union_points);
  CHECK(BoundaryPhysicalTraceCache::LocalUnionGroupCount() ==
        local_before + expected_local_groups);
  const long long local_operators =
      BoundaryPhysicalTraceCache::LocalUnionOperatorCount() - local_operators_before;
  CHECK(local_operators <= expected_local_groups);
  CHECK((expected_local_groups == 0 || local_operators > 0));
  CHECK(BoundaryPhysicalTraceCache::GhostUnionGroupCount() ==
        ghost_before + expected_ghost_groups);
  CHECK(BoundaryPhysicalTraceCache::RequestCountBefore() ==
        requests_before + expected_requests_before);
  CHECK(BoundaryPhysicalTraceCache::RequestCountAfter() ==
        requests_after + expected_requests_after);
  CHECK(BoundaryPhysicalTraceCache::UnionRouteCount() == routes_before + expected_routes);
  CHECK(BoundaryPhysicalTraceCache::RoutingOperatorCount() ==
        routing_operators_before +
            static_cast<long long>(expected_routing_patterns.size()));
  if (expect_ghost_union)
  {
    long long ghost_deltas[3] = {
        BoundaryPhysicalTraceCache::GhostUnionGroupCount() - ghost_before,
        BoundaryPhysicalTraceCache::RequestCountBefore() - requests_before,
        BoundaryPhysicalTraceCache::RequestCountAfter() - requests_after};
    Mpi::GlobalSum(3, ghost_deltas, comm);
    CHECK(ghost_deltas[0] > 0);
    CHECK(ghost_deltas[1] > ghost_deltas[2]);
    CHECK(ghost_deltas[2] > 0);
  }
}

}  // namespace

TEST_CASE("DefaultIntegrationOrder Periodic Trace", "[surfacefunctional][Serial]")
{
  mfem::Mesh orig_mesh = mfem::Mesh::MakeCartesian3D(3, 3, 3, mfem::Element::HEXAHEDRON);
  orig_mesh.EnsureNodes();
  std::vector<mfem::Vector> translations = {mfem::Vector({1.0, 0.0, 0.0}),
                                            mfem::Vector({0.0, 1.0, 0.0})};
  mfem::Mesh mesh = mfem::Mesh::MakePeriodic(
      orig_mesh, orig_mesh.CreatePeriodicVertexMapping(translations));

  REQUIRE(mesh.GetNBE() > 0);
  REQUIRE(mesh.GetNodes());
  const auto geom = mesh.GetBdrElementGeometry(0);
  const auto *fec = mesh.GetNodalFESpace()->FEColl();
  CHECK(fec->FiniteElementForGeometry(geom) == nullptr);
  REQUIRE(fec->TraceFiniteElementForGeometry(geom));

  fem::DefaultIntegrationOrder::p_trial = 1;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;
  const auto *T = mesh.GetBdrElementTransformation(0);
  REQUIRE(T);
  CHECK(fem::DefaultIntegrationOrder::Get(mesh, geom) ==
        fem::DefaultIntegrationOrder::Get(*T));
}

TEST_CASE("SurfaceFunctional Area", "[surfacefunctional][Serial][Parallel][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  fem::DefaultIntegrationOrder::p_trial = 1;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;
  auto file =
      GENERATE(std::string("fichera-tet.mesh"), std::string("fichera-hex.mesh"),
               std::string("fichera-mixed-p2.mesh"), std::string("tinyzoo-3d.mesh"));
  CAPTURE(file);
  auto mesh = LoadTestMesh(comm, file);
  auto &pmesh = mesh->Get();

  const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;

  SECTION("All boundary attributes")
  {
    mfem::Array<int> marker(bdr_attr_max);
    marker = 1;
    SurfaceFunctional area(SurfaceFunctional::Kind::AREA, *mesh, marker);
    CHECK(area.Eval() == Catch::Approx(RefSurfaceArea(pmesh, marker)).epsilon(1.0e-12));
    if (file != "fichera-mixed-p2.mesh")
    {
      // For meshes with straight-sided boundary elements, also cross-check against the
      // independent implementation in mesh::GetSurfaceArea (exact for any quadrature
      // order, so quadrature rule differences don't matter).
      CHECK(area.Eval() ==
            Catch::Approx(mesh::GetSurfaceArea(pmesh, marker)).epsilon(1.0e-12));
    }
  }

  SECTION("Single boundary attribute")
  {
    for (auto attr : pmesh.bdr_attributes)
    {
      CAPTURE(attr);
      mfem::Array<int> marker(bdr_attr_max);
      marker = 0;
      marker[attr - 1] = 1;
      SurfaceFunctional area(SurfaceFunctional::Kind::AREA, *mesh, marker);
      CHECK(area.Eval() == Catch::Approx(RefSurfaceArea(pmesh, marker)).epsilon(1.0e-12));
    }
  }
}

TEST_CASE("SurfaceFunctional HCurl Norm (fixed mapped face rules)",
          "[surfacefunctional][Serial][Parallel][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto file =
      GENERATE(std::string("fichera-tet.mesh"), std::string("fichera-hex.mesh"),
               std::string("fichera-mixed-p2.mesh"), std::string("tinyzoo-3d.mesh"));
  auto order = GENERATE(1, 2);
  CAPTURE(file, order);
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;
  auto mesh = LoadTestMesh(comm, file);
  auto &pmesh = mesh->Get();
  if (file == "tinyzoo-3d.mesh")
  {
    // The mixed zoo has boundary faces on both the prism and pyramid. The all-boundary
    // section below therefore exercises their fixed mapped volume rules, alongside the
    // existing tet, hex, and mixed-mesh cases in this generator.
    bool has_prism = false, has_pyramid = false;
    for (int i = 0; i < pmesh.GetNBE(); i++)
    {
      int elem, face_info;
      pmesh.GetBdrElementAdjacentElement(i, elem, face_info);
      const auto geom = pmesh.GetElementGeometry(elem);
      has_prism = has_prism || geom == mfem::Geometry::PRISM;
      has_pyramid = has_pyramid || geom == mfem::Geometry::PYRAMID;
    }
    Mpi::GlobalOr(1, &has_prism, comm);
    Mpi::GlobalOr(1, &has_pyramid, comm);
    REQUIRE(has_prism);
    REQUIRE(has_pyramid);
  }

  mfem::ND_FECollection nd_fec(order, pmesh.Dimension());
  FiniteElementSpace nd_fespace(*mesh, &nd_fec);

  // Project a non-trivial smooth vector field onto the ND space. The reference value is
  // computed from the same projected grid function, so the comparison is exact up to
  // quadrature evaluation differences (not projection error).
  mfem::ParGridFunction E(&nd_fespace.Get());
  mfem::VectorFunctionCoefficient f(3,
                                    [](const mfem::Vector &x, mfem::Vector &v)
                                    {
                                      v(0) = std::sin(x(1)) + x(2) * x(2);
                                      v(1) = std::cos(x(2)) + x(0);
                                      v(2) = x(0) * x(1) + 1.0;
                                    });
  E.ProjectCoefficient(f);

  const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;

  SECTION("All boundary attributes")
  {
    mfem::Array<int> marker(bdr_attr_max);
    marker = 1;
    SurfaceFunctional norm2(SurfaceFunctional::Kind::HCURL_NORM2, *mesh, marker,
                            &nd_fespace.Get());
    const double ref = RefSurfaceHCurlNorm2(pmesh, E, marker);
    CHECK(norm2.Eval(&E) == Catch::Approx(ref).epsilon(1.0e-10));
  }

  SECTION("Single boundary attribute")
  {
    for (auto attr : pmesh.bdr_attributes)
    {
      CAPTURE(attr);
      mfem::Array<int> marker(bdr_attr_max);
      marker = 0;
      marker[attr - 1] = 1;
      SurfaceFunctional norm2(SurfaceFunctional::Kind::HCURL_NORM2, *mesh, marker,
                              &nd_fespace.Get());
      const double ref = RefSurfaceHCurlNorm2(pmesh, E, marker);
      CHECK(norm2.Eval(&E) == Catch::Approx(ref).epsilon(1.0e-10));
    }
  }
}

TEST_CASE("SurfaceFunctional Interface Dielectric", "[surfacefunctional][Serial][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto elem_type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  auto order = GENERATE(1, 2);
  auto complex = GENERATE(false, true);
  CAPTURE(elem_type, order, complex);
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  auto mesh = MakeInterfaceMesh(comm, elem_type);
  auto &pmesh = mesh->Get();

  // Materials: vacuum for attribute 1, dielectric for attribute 2.
  config::MaterialData vacuum, dielectric;
  vacuum.attributes = {1};
  dielectric.attributes = {2};
  dielectric.epsilon_r.s[0] = 11.7;
  dielectric.epsilon_r.s[1] = 11.7;
  dielectric.epsilon_r.s[2] = 11.7;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({vacuum, dielectric}, periodic, ProblemType::DRIVEN, *mesh);

  mfem::ND_FECollection nd_fec(order, pmesh.Dimension());
  FiniteElementSpace nd_fespace(*mesh, &nd_fec);

  // Project non-trivial smooth vector fields onto the (complex-valued) ND space.
  GridFunction E(nd_fespace, complex);
  mfem::VectorFunctionCoefficient fr(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = std::sin(x(1)) + x(2) * x(2);
                                       v(1) = std::cos(x(2)) + x(0);
                                       v(2) = x(0) * x(1) + 1.0;
                                     });
  E.Real().ProjectCoefficient(fr);
  if (complex)
  {
    mfem::VectorFunctionCoefficient fi(3,
                                       [](const mfem::Vector &x, mfem::Vector &v)
                                       {
                                         v(0) = x(1) * x(2) - 0.5;
                                         v(1) = std::sin(x(0)) - x(2);
                                         v(2) = std::cos(x(1)) + x(0) * x(0);
                                       });
    E.Imag().ProjectCoefficient(fi);
  }

  const double t_i = 2.0e-3, epsilon_i = 10.0;
  const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;

  auto TestType = [&](InterfaceDielectric type, const mfem::Array<int> &marker)
  {
    SurfaceFunctional epr(*mesh, marker, nd_fespace, mat_op, type, t_i, epsilon_i);
    auto MakeLegacy = [&]() -> std::unique_ptr<mfem::Coefficient>
    {
      switch (type)
      {
        case InterfaceDielectric::DEFAULT:
          return std::make_unique<
              InterfaceDielectricCoefficient<InterfaceDielectric::DEFAULT>>(E, mat_op, t_i,
                                                                            epsilon_i);
        case InterfaceDielectric::MA:
          return std::make_unique<InterfaceDielectricCoefficient<InterfaceDielectric::MA>>(
              E, mat_op, t_i, epsilon_i);
        case InterfaceDielectric::MS:
          return std::make_unique<InterfaceDielectricCoefficient<InterfaceDielectric::MS>>(
              E, mat_op, t_i, epsilon_i);
        case InterfaceDielectric::SA:
          return std::make_unique<InterfaceDielectricCoefficient<InterfaceDielectric::SA>>(
              E, mat_op, t_i, epsilon_i);
      }
      return {};
    };
    auto legacy = MakeLegacy();
    const double ref = RefSurfaceCoefficientIntegral(pmesh, *legacy, marker);
    const double val = epr.Eval(E);
    CAPTURE(ref, val);
    if (std::abs(ref) > 0.0)
    {
      CHECK(val == Catch::Approx(ref).epsilon(1.0e-10));
    }
    else
    {
      CHECK(std::abs(val) < 1.0e-12);
    }
  };

  SECTION("Exterior boundary")
  {
    // Boundary attribute 1 is an exterior boundary (single-sided evaluation). The
    // bottom boundary (z = 0) neighbors vacuum, the top (z = 1) neighbors dielectric.
    for (auto type : {InterfaceDielectric::DEFAULT, InterfaceDielectric::MA,
                      InterfaceDielectric::MS, InterfaceDielectric::SA})
    {
      CAPTURE(static_cast<int>(type));
      for (int attr : {1, 6})  // z = 0 (vacuum side), z = 1 (dielectric side)
      {
        CAPTURE(attr);
        mfem::Array<int> marker(bdr_attr_max);
        marker = 0;
        marker[attr - 1] = 1;
        TestType(type, marker);
      }
    }
  }

  SECTION("Interior material interface")
  {
    // Boundary attribute 7 is the interior vacuum-dielectric interface (two-sided,
    // side selection by material light speed for MA/MS/SA, averaging for DEFAULT).
    mfem::Array<int> marker(bdr_attr_max);
    marker = 0;
    marker[7 - 1] = 1;
    for (auto type : {InterfaceDielectric::DEFAULT, InterfaceDielectric::MA,
                      InterfaceDielectric::MS, InterfaceDielectric::SA})
    {
      CAPTURE(static_cast<int>(type));
      TestType(type, marker);
    }
  }
}

TEST_CASE("SurfaceFunctional Interface Dielectric 2D", "[surfacefunctional][Serial][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto elem_type = GENERATE(mfem::Element::TRIANGLE, mfem::Element::QUADRILATERAL);
  auto order = GENERATE(1, 2);
  auto complex = GENERATE(false, true);
  CAPTURE(elem_type, order, complex);
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  auto mesh = MakeInterfaceMesh2D(comm, elem_type);
  auto &pmesh = mesh->Get();

  // Materials: vacuum for attribute 1, dielectric for attribute 2. Use diagonal
  // anisotropy to exercise the 2x2 coefficient unpacking path in the libCEED kernels.
  config::MaterialData vacuum, dielectric;
  vacuum.attributes = {1};
  dielectric.attributes = {2};
  dielectric.epsilon_r.s[0] = 11.7;
  dielectric.epsilon_r.s[1] = 3.2;
  dielectric.epsilon_r.s[2] = 5.1;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({vacuum, dielectric}, periodic, ProblemType::DRIVEN, *mesh);

  mfem::ND_FECollection nd_fec(order, pmesh.Dimension());
  FiniteElementSpace nd_fespace(*mesh, &nd_fec);

  GridFunction E(nd_fespace, complex);
  mfem::VectorFunctionCoefficient fr(2,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = std::sin(1.7 * x(1)) + 0.25 * x(0);
                                       v(1) = std::cos(0.3 + x(0)) + x(1) * x(1);
                                     });
  E.Real().ProjectCoefficient(fr);
  if (complex)
  {
    mfem::VectorFunctionCoefficient fi(2,
                                       [](const mfem::Vector &x, mfem::Vector &v)
                                       {
                                         v(0) = x(0) * x(1) - 0.2;
                                         v(1) = std::sin(x(0) + 0.5 * x(1));
                                       });
    E.Imag().ProjectCoefficient(fi);
  }

  const double t_i = 2.0e-3, epsilon_i = 10.0;
  const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;

  auto TestType = [&](InterfaceDielectric type, const mfem::Array<int> &marker)
  {
    SurfaceFunctional epr(*mesh, marker, nd_fespace, mat_op, type, t_i, epsilon_i);
    REQUIRE(epr.IsValid());
    auto MakeLegacy = [&]() -> std::unique_ptr<mfem::Coefficient>
    {
      switch (type)
      {
        case InterfaceDielectric::DEFAULT:
          return std::make_unique<
              InterfaceDielectricCoefficient<InterfaceDielectric::DEFAULT>>(E, mat_op, t_i,
                                                                            epsilon_i);
        case InterfaceDielectric::MA:
          return std::make_unique<InterfaceDielectricCoefficient<InterfaceDielectric::MA>>(
              E, mat_op, t_i, epsilon_i);
        case InterfaceDielectric::MS:
          return std::make_unique<InterfaceDielectricCoefficient<InterfaceDielectric::MS>>(
              E, mat_op, t_i, epsilon_i);
        case InterfaceDielectric::SA:
          return std::make_unique<InterfaceDielectricCoefficient<InterfaceDielectric::SA>>(
              E, mat_op, t_i, epsilon_i);
      }
      return {};
    };
    auto legacy = MakeLegacy();
    const double ref = RefSurfaceCoefficientIntegral(pmesh, *legacy, marker);
    const double val = epr.Eval(E);
    CAPTURE(ref, val);
    if (std::abs(ref) > 0.0)
    {
      CHECK(val == Catch::Approx(ref).epsilon(1.0e-10));
    }
    else
    {
      CHECK(std::abs(val) < 1.0e-12);
    }
  };

  SECTION("Exterior boundary")
  {
    for (auto type : {InterfaceDielectric::DEFAULT, InterfaceDielectric::MA,
                      InterfaceDielectric::MS, InterfaceDielectric::SA})
    {
      CAPTURE(static_cast<int>(type));
      for (int attr : {1, 3})  // bottom (vacuum side), top (dielectric side)
      {
        CAPTURE(attr);
        REQUIRE(pmesh.bdr_attributes.Find(attr) >= 0);
        mfem::Array<int> marker(bdr_attr_max);
        marker = 0;
        marker[attr - 1] = 1;
        TestType(type, marker);
      }
    }
  }

  SECTION("Interior material interface")
  {
    mfem::Array<int> marker(bdr_attr_max);
    marker = 0;
    marker[7 - 1] = 1;
    for (auto type : {InterfaceDielectric::DEFAULT, InterfaceDielectric::MA,
                      InterfaceDielectric::MS, InterfaceDielectric::SA})
    {
      CAPTURE(static_cast<int>(type));
      TestType(type, marker);
    }
  }
}

TEST_CASE("SurfaceFunctional 2D nonconformal mapped subfaces",
          "[surfacefunctional][Serial][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const auto elem_type = GENERATE(mfem::Element::TRIANGLE, mfem::Element::QUADRILATERAL);
  const int refinement_depth = GENERATE(1, 2);
  CAPTURE(elem_type, refinement_depth);
  fem::DefaultIntegrationOrder::p_trial = 1;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  auto mesh = MakeNCInterfaceMesh2D(comm, elem_type, refinement_depth);
  auto &pmesh = mesh->Get();
  REQUIRE(pmesh.Nonconforming());
  REQUIRE(pmesh.bdr_attributes.Find(7) >= 0);

  config::MaterialData vacuum, dielectric;
  vacuum.attributes = {1};
  dielectric.attributes = {2};
  dielectric.epsilon_r.s[0] = 11.7;
  dielectric.epsilon_r.s[1] = 3.2;
  dielectric.epsilon_r.s[2] = 5.1;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({vacuum, dielectric}, periodic, ProblemType::DRIVEN, *mesh);

  mfem::ND_FECollection nd_fec(1, 2);
  FiniteElementSpace nd_fespace(*mesh, &nd_fec);
  GridFunction E(nd_fespace);
  mfem::VectorFunctionCoefficient field(2,
                                        [](const mfem::Vector &x, mfem::Vector &v)
                                        {
                                          v(0) = std::sin(1.7 * x(1)) + 0.25 * x(0);
                                          v(1) = std::cos(0.3 + x(0)) + x(1) * x(1);
                                        });
  E.Real().ProjectCoefficient(field);

  mfem::Array<int> marker(pmesh.bdr_attributes.Max());
  marker = 0;
  marker[7 - 1] = 1;
  constexpr double t_i = 2.0e-3, epsilon_i = 10.0;
  SurfaceFunctional epr(*mesh, marker, nd_fespace, mat_op, InterfaceDielectric::DEFAULT,
                        t_i, epsilon_i);
  REQUIRE(epr.IsValid());
  InterfaceDielectricCoefficient<InterfaceDielectric::DEFAULT> legacy(E, mat_op, t_i,
                                                                      epsilon_i);
  const double ref = RefSurfaceCoefficientIntegral(pmesh, legacy, marker);
  CHECK(epr.Eval(E) == Catch::Approx(ref).epsilon(1.0e-10));
}

TEST_CASE_METHOD(test::SharedTempDir,
                 "PostOperator boundary mode 2D interface dielectric matches legacy",
                 "[postoperator][surfacefunctional][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  REQUIRE(Mpi::Size(comm) == 1);
  constexpr int order = 2;

  IoData iodata(Units(1.0, 1.0));
  iodata.problem.type = ProblemType::BOUNDARYMODE;
  iodata.problem.verbose = 0;
  iodata.problem.output = temp_dir.string();
  iodata.problem.output_formats.paraview = false;
  iodata.problem.output_formats.gridfunction = false;
  iodata.model.L0 = 1.0;
  iodata.model.Lc = 1.0;
  iodata.solver.order = order;
  iodata.solver.boundary_mode.freq = 1.0;
  iodata.solver.boundary_mode.n = 1;
  iodata.solver.boundary_mode.n_post = 1;
  iodata.solver.linear.tol = 1.0e-8;
  iodata.solver.linear.max_it = 50;

  config::MaterialData vacuum, dielectric;
  vacuum.attributes = {1};
  dielectric.attributes = {2};
  dielectric.epsilon_r.s[0] = 4.0;
  dielectric.epsilon_r.s[1] = 6.0;
  dielectric.epsilon_r.s[2] = 8.0;
  iodata.domains.materials = {vacuum, dielectric};
  iodata.boundaries.pec.attributes = {1, 2, 3, 4};

  config::InterfaceDielectricData interface;
  interface.attributes = {7};
  interface.type = InterfaceDielectric::SA;
  interface.t = 2.0e-3;
  interface.epsilon_r = 10.0;
  interface.tandelta = 3.0e-4;
  iodata.boundaries.postpro.dielectric.emplace(1, interface);

  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(MakeInterfaceMesh2D(comm, mfem::Element::TRIANGLE));
  auto &pmesh = mesh.front()->Get();
  MaterialOperator mat_op(iodata, *mesh.front());
  BoundaryModeOperator fem_op(iodata, mesh, mat_op);
  TestBoundaryModePostOperator post_op(iodata, fem_op);

  GridFunction E(fem_op.GetNDSpace(), true), En(fem_op.GetH1Space(), true);
  mfem::VectorFunctionCoefficient etr(2,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::sin(x(1)) + 0.3 * x(0);
                                        v(1) = std::cos(x(0)) + x(0) * x(1);
                                      });
  mfem::VectorFunctionCoefficient eti(2,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(0) - 0.2 * x(1);
                                        v(1) = std::sin(x(0) + x(1));
                                      });
  mfem::FunctionCoefficient enr([](const mfem::Vector &x) { return 0.5 + x(0) * x(1); });
  mfem::FunctionCoefficient eni([](const mfem::Vector &x)
                                { return std::cos(x(0)) - 0.1 * x(1); });
  E.Real().ProjectCoefficient(etr);
  E.Imag().ProjectCoefficient(eti);
  En.Real().ProjectCoefficient(enr);
  En.Imag().ProjectCoefficient(eni);

  ComplexVector et(fem_op.GetNDTrueVSize()), en(fem_op.GetH1TrueVSize());
  E.Real().GetTrueDofs(et.Real());
  E.Imag().GetTrueDofs(et.Imag());
  En.Real().GetTrueDofs(en.Real());
  En.Imag().GetTrueDofs(en.Imag());

  // Public PostOperator entry point under test: BoundaryMode computes its derived B
  // fields internally, then measures configured interface dielectric participation.
  const double omega = 2.0;
  const std::complex<double> kn(0.7, 0.05);
  post_op.MeasureAndPrintAll(0, et, en, kn, omega, 0.0, 0.0, 1);

  mfem::Array<int> marker(pmesh.bdr_attributes.Max());
  marker = 0;
  marker[7 - 1] = 1;
  InterfaceDielectricCoefficient<InterfaceDielectric::SA> legacy(E, mat_op, interface.t,
                                                                 interface.epsilon_r);
  const double interface_energy = RefSurfaceCoefficientIntegral(pmesh, legacy, marker);
  DomainPostOperator domain_post(iodata.domains.postpro, mat_op, fem_op.GetNDSpace(),
                                 fem_op.GetCurlSpace());
  const double domain_energy = domain_post.GetElectricFieldEnergy(E);
  const double participation = interface_energy / domain_energy;

  const auto &cache = post_op.Cache();
  REQUIRE(cache.interface_eps_i.size() == 1);
  const auto &measured = cache.interface_eps_i.front();
  CHECK(measured.idx == 1);
  CHECK(measured.energy == Catch::Approx(interface_energy).epsilon(1.0e-10));
  CHECK(measured.energy_participation == Catch::Approx(participation).epsilon(1.0e-10));
  CHECK(measured.quality_factor ==
        Catch::Approx(1.0 / (participation * interface.tandelta)).epsilon(1.0e-10));
}

TEST_CASE("SurfaceFunctional Surface Flux", "[surfacefunctional][Serial][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto elem_type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  auto order = GENERATE(1, 2);
  auto complex = GENERATE(false, true);
  CAPTURE(elem_type, order, complex);
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  auto mesh = MakeInterfaceMesh(comm, elem_type);
  auto &pmesh = mesh->Get();

  config::MaterialData vacuum, dielectric;
  vacuum.attributes = {1};
  dielectric.attributes = {2};
  dielectric.epsilon_r.s[0] = 11.7;
  dielectric.epsilon_r.s[1] = 11.7;
  dielectric.epsilon_r.s[2] = 11.7;
  dielectric.mu_r.s[0] = 1.4;
  dielectric.mu_r.s[1] = 1.4;
  dielectric.mu_r.s[2] = 1.4;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({vacuum, dielectric}, periodic, ProblemType::DRIVEN, *mesh);

  mfem::ND_FECollection nd_fec(order, pmesh.Dimension());
  mfem::RT_FECollection rt_fec(order - 1, pmesh.Dimension());
  FiniteElementSpace nd_fespace(*mesh, &nd_fec), rt_fespace(*mesh, &rt_fec);

  // Project non-trivial smooth vector fields onto the (complex-valued) ND and RT
  // spaces.
  GridFunction E(nd_fespace, complex), B(rt_fespace, complex);
  mfem::VectorFunctionCoefficient fer(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::sin(x(1)) + x(2) * x(2);
                                        v(1) = std::cos(x(2)) + x(0);
                                        v(2) = x(0) * x(1) + 1.0;
                                      });
  mfem::VectorFunctionCoefficient fbr(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) - 0.3 * x(2);
                                        v(1) = std::sin(x(2)) + 0.5;
                                        v(2) = std::cos(x(0)) - x(1) * x(2);
                                      });
  E.Real().ProjectCoefficient(fer);
  B.Real().ProjectCoefficient(fbr);
  if (complex)
  {
    mfem::VectorFunctionCoefficient fei(3,
                                        [](const mfem::Vector &x, mfem::Vector &v)
                                        {
                                          v(0) = x(1) * x(2) - 0.5;
                                          v(1) = std::sin(x(0)) - x(2);
                                          v(2) = std::cos(x(1)) + x(0) * x(0);
                                        });
    mfem::VectorFunctionCoefficient fbi(3,
                                        [](const mfem::Vector &x, mfem::Vector &v)
                                        {
                                          v(0) = std::cos(x(2)) - 0.2;
                                          v(1) = x(0) * x(2) + 0.1;
                                          v(2) = std::sin(x(1)) - x(0);
                                        });
    E.Imag().ProjectCoefficient(fei);
    B.Imag().ProjectCoefficient(fbi);
  }

  mfem::Vector x0(3);
  x0(0) = 0.4;
  x0(1) = 0.6;
  x0(2) = -0.2;  // Off-center so the orientation sign flip logic is exercised
  const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;

  auto TestType = [&](SurfaceFlux type, bool two_sided, const mfem::Array<int> &marker)
  {
    const GridFunction *E_use =
        (type == SurfaceFlux::ELECTRIC || type == SurfaceFlux::POWER) ? &E : nullptr;
    const GridFunction *B_use =
        (type == SurfaceFlux::MAGNETIC || type == SurfaceFlux::POWER) ? &B : nullptr;
    SurfaceFunctional flux(*mesh, marker, E_use ? &nd_fespace.Get() : nullptr,
                           B_use ? &rt_fespace.Get() : nullptr, mat_op, type, two_sided,
                           x0);

    // Legacy reference following SurfacePostOperator::GetSurfaceFlux.
    auto MakeLegacy =
        [&](const mfem::ParGridFunction *Er,
            const mfem::ParGridFunction *Br) -> std::unique_ptr<mfem::Coefficient>
    {
      switch (type)
      {
        case SurfaceFlux::ELECTRIC:
          return std::make_unique<BdrSurfaceFluxCoefficient<SurfaceFlux::ELECTRIC>>(
              Er, Br, mat_op, two_sided, x0);
        case SurfaceFlux::MAGNETIC:
          return std::make_unique<BdrSurfaceFluxCoefficient<SurfaceFlux::MAGNETIC>>(
              Er, Br, mat_op, two_sided, x0);
        case SurfaceFlux::POWER:
          return std::make_unique<BdrSurfaceFluxCoefficient<SurfaceFlux::POWER>>(
              Er, Br, mat_op, two_sided, x0);
      }
      return {};
    };
    auto legacy_re = MakeLegacy(E_use ? &E.Real() : nullptr, B_use ? &B.Real() : nullptr);
    std::complex<double> ref(RefSurfaceCoefficientIntegral(pmesh, *legacy_re, marker), 0.0);
    if (complex)
    {
      auto legacy_im = MakeLegacy(E_use ? &E.Imag() : nullptr, B_use ? &B.Imag() : nullptr);
      const double ref_im = RefSurfaceCoefficientIntegral(pmesh, *legacy_im, marker);
      if (type == SurfaceFlux::POWER)
      {
        ref += ref_im;
      }
      else
      {
        ref.imag(ref_im);
      }
    }

    const std::complex<double> val = flux.EvalFlux(E_use, B_use);
    CAPTURE(ref.real(), ref.imag(), val.real(), val.imag());
    auto CheckPart = [](double v, double r)
    {
      if (std::abs(r) > 1.0e-12)
      {
        CHECK(v == Catch::Approx(r).epsilon(1.0e-10));
      }
      else
      {
        CHECK(std::abs(v) < 1.0e-10);
      }
    };
    CheckPart(val.real(), ref.real());
    CheckPart(val.imag(), ref.imag());
  };

  SECTION("Exterior boundary")
  {
    for (auto type : {SurfaceFlux::ELECTRIC, SurfaceFlux::MAGNETIC, SurfaceFlux::POWER})
    {
      CAPTURE(static_cast<int>(type));
      for (int attr : {1, 6})  // z = 0 (vacuum side), z = 1 (dielectric side)
      {
        CAPTURE(attr);
        mfem::Array<int> marker(bdr_attr_max);
        marker = 0;
        marker[attr - 1] = 1;
        TestType(type, false, marker);
      }
    }
  }

  SECTION("Interior material interface")
  {
    mfem::Array<int> marker(bdr_attr_max);
    marker = 0;
    marker[7 - 1] = 1;
    for (auto type : {SurfaceFlux::ELECTRIC, SurfaceFlux::MAGNETIC, SurfaceFlux::POWER})
    {
      CAPTURE(static_cast<int>(type));
      for (bool two_sided : {false, true})
      {
        CAPTURE(two_sided);
        TestType(type, two_sided, marker);
      }
    }
  }
}

TEST_CASE("SurfaceFunctional Complex Power", "[surfacefunctional][Serial][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto elem_type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  auto order = GENERATE(1, 2);
  auto complex = GENERATE(false, true);
  CAPTURE(elem_type, order, complex);
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  auto mesh = MakeInterfaceMesh(comm, elem_type);
  auto &pmesh = mesh->Get();

  config::MaterialData vacuum, dielectric;
  vacuum.attributes = {1};
  dielectric.attributes = {2};
  dielectric.epsilon_r.s[0] = 11.7;
  dielectric.epsilon_r.s[1] = 11.7;
  dielectric.epsilon_r.s[2] = 11.7;
  dielectric.mu_r.s[0] = 1.4;
  dielectric.mu_r.s[1] = 1.4;
  dielectric.mu_r.s[2] = 1.4;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({vacuum, dielectric}, periodic, ProblemType::DRIVEN, *mesh);

  mfem::ND_FECollection nd_fec(order, pmesh.Dimension());
  mfem::RT_FECollection rt_fec(order - 1, pmesh.Dimension());
  FiniteElementSpace nd_fespace(*mesh, &nd_fec), rt_fespace(*mesh, &rt_fec);

  GridFunction E(nd_fespace, complex), B(rt_fespace, complex);
  mfem::VectorFunctionCoefficient fer(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::sin(x(1)) + x(2) * x(2);
                                        v(1) = std::cos(x(2)) + x(0);
                                        v(2) = x(0) * x(1) + 1.0;
                                      });
  mfem::VectorFunctionCoefficient fbr(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) - 0.3 * x(2);
                                        v(1) = std::sin(x(2)) + 0.5;
                                        v(2) = std::cos(x(0)) - x(1) * x(2);
                                      });
  E.Real().ProjectCoefficient(fer);
  B.Real().ProjectCoefficient(fbr);
  if (complex)
  {
    mfem::VectorFunctionCoefficient fei(3,
                                        [](const mfem::Vector &x, mfem::Vector &v)
                                        {
                                          v(0) = x(1) * x(2) - 0.5;
                                          v(1) = std::sin(x(0)) - x(2);
                                          v(2) = std::cos(x(1)) + x(0) * x(0);
                                        });
    mfem::VectorFunctionCoefficient fbi(3,
                                        [](const mfem::Vector &x, mfem::Vector &v)
                                        {
                                          v(0) = std::cos(x(2)) - 0.2;
                                          v(1) = x(0) * x(2) + 0.1;
                                          v(2) = std::sin(x(1)) - x(0);
                                        });
    E.Imag().ProjectCoefficient(fei);
    B.Imag().ProjectCoefficient(fbi);
  }

  // Legacy reference: replicate LumpedPortData::GetPower exactly (linear form from
  // BdrSurfaceCurrentVectorCoefficient dotted with the E field expansion).
  auto LegacyPower = [&](int attr, const mfem::Array<int> &attr_marker)
  {
    mfem::Array<int> attr_list(1);
    attr_list[0] = attr;
    std::complex<double> dot;
    {
      RestrictedVectorCoefficient<BdrSurfaceCurrentVectorCoefficient> fbr_c(
          attr_list, B.Real(), mat_op);
      mfem::LinearForm pr(&nd_fespace.Get());
      pr.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fbr_c),
                               const_cast<mfem::Array<int> &>(attr_marker));
      pr.UseFastAssembly(false);
      pr.UseDevice(false);
      pr.Assemble();
      pr.UseDevice(true);
      dot =
          -(pr * E.Real()) + (complex ? std::complex<double>(0.0, -(pr * E.Imag())) : 0.0);
    }
    if (complex)
    {
      RestrictedVectorCoefficient<BdrSurfaceCurrentVectorCoefficient> fbi_c(
          attr_list, B.Imag(), mat_op);
      mfem::LinearForm pi(&nd_fespace.Get());
      pi.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fbi_c),
                               const_cast<mfem::Array<int> &>(attr_marker));
      pi.UseFastAssembly(false);
      pi.UseDevice(false);
      pi.Assemble();
      pi.UseDevice(true);
      dot += std::complex<double>(-(pi * E.Imag()), pi * E.Real());
    }
    Mpi::GlobalSum(1, &dot, comm);
    return dot;
  };

  const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
  mfem::Vector x0(3);
  x0 = 0.0;

  auto CheckPart = [](double v, double r)
  {
    if (std::abs(r) > 1.0e-12)
    {
      CHECK(v == Catch::Approx(r).epsilon(1.0e-10));
    }
    else
    {
      CHECK(std::abs(v) < 1.0e-10);
    }
  };

  for (int attr : {1, 6, 7})  // Exterior (vacuum side), exterior (dielectric), interior
  {
    CAPTURE(attr);
    mfem::Array<int> marker(bdr_attr_max);
    marker = 0;
    marker[attr - 1] = 1;
    SurfaceFunctional power(*mesh, marker, &nd_fespace.Get(), &rt_fespace.Get(), mat_op,
                            SurfaceFlux::POWER, /*two_sided*/ true, x0);
    const std::complex<double> ref = LegacyPower(attr, marker);
    const std::complex<double> val = power.EvalComplexPower(E, B);
    CAPTURE(ref.real(), ref.imag(), val.real(), val.imag());
    CheckPart(val.real(), ref.real());
    CheckPart(val.imag(), ref.imag());
  }

  SECTION("Batched by boundary attribute")
  {
    mfem::Array<int> marker(bdr_attr_max), attr_to_bin(bdr_attr_max);
    marker = 0;
    attr_to_bin = -1;
    std::vector<int> attrs = {1, 6, 7};
    for (std::size_t i = 0; i < attrs.size(); i++)
    {
      marker[attrs[i] - 1] = 1;
      attr_to_bin[attrs[i] - 1] = static_cast<int>(i);
    }
    SurfaceFunctional power(*mesh, marker, &nd_fespace.Get(), &rt_fespace.Get(), mat_op,
                            SurfaceFlux::POWER, /*two_sided*/ true, x0);
    const auto values = power.EvalComplexPowerByAttribute(E, B, attr_to_bin,
                                                          static_cast<int>(attrs.size()));
    REQUIRE(values.size() == attrs.size());
    for (std::size_t i = 0; i < attrs.size(); i++)
    {
      mfem::Array<int> single_marker(bdr_attr_max);
      single_marker = 0;
      single_marker[attrs[i] - 1] = 1;
      const std::complex<double> ref = LegacyPower(attrs[i], single_marker);
      CAPTURE(attrs[i], ref.real(), ref.imag(), values[i].real(), values[i].imag());
      CheckPart(values[i].real(), ref.real());
      CheckPart(values[i].imag(), ref.imag());
    }
  }
}

// Benchmark comparing the legacy mfem::Coefficient-based measurement paths against the
// libCEED surface functionals. Hidden by default; run explicitly with:
//   ./palace-unit-tests "[surfacefunctional-bench]" [--device cuda]
// Mesh size can be controlled with the PALACE_BENCH_N environment variable (default
// 20, i.e. a 20x20x20 hex mesh).
TEST_CASE("SurfaceFunctional Benchmark", "[.][surfacefunctional-bench][Serial][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const int order = 2;
  const int N =
      std::getenv("PALACE_BENCH_N") ? std::atoi(std::getenv("PALACE_BENCH_N")) : 20;
  const int n_reps = 10;
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  // Two-material N x N x N hex mesh with an interior interface (attribute 7).
  auto mesh = [&]()
  {
    mfem::Mesh smesh = mfem::Mesh::MakeCartesian3D(N, N, N, mfem::Element::HEXAHEDRON);
    for (int e = 0; e < smesh.GetNE(); e++)
    {
      mfem::Vector center(3);
      smesh.GetElementCenter(e, center);
      smesh.SetAttribute(e, (center(2) < 0.5) ? 1 : 2);
    }
    for (int f = 0; f < smesh.GetNumFaces(); f++)
    {
      int e1, e2;
      smesh.GetFaceElements(f, &e1, &e2);
      if (e1 >= 0 && e2 >= 0 && smesh.GetAttribute(e1) != smesh.GetAttribute(e2))
      {
        auto *face_elem = smesh.GetFace(f)->Duplicate(&smesh);
        face_elem->SetAttribute(7);
        smesh.AddBdrElement(face_elem);
      }
    }
    smesh.FinalizeTopology();
    smesh.Finalize();
    smesh.SetAttributes();
    smesh.EnsureNodes();
    auto pmesh = std::make_unique<mfem::ParMesh>(comm, smesh);
    return std::make_unique<Mesh>(std::move(pmesh));
  }();
  auto &pmesh = mesh->Get();

  config::MaterialData vacuum, dielectric;
  vacuum.attributes = {1};
  dielectric.attributes = {2};
  dielectric.epsilon_r.s[0] = 11.7;
  dielectric.epsilon_r.s[1] = 11.7;
  dielectric.epsilon_r.s[2] = 11.7;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({vacuum, dielectric}, periodic, ProblemType::DRIVEN, *mesh);

  mfem::ND_FECollection nd_fec(order, 3);
  mfem::RT_FECollection rt_fec(order - 1, 3);
  mfem::H1_FECollection h1_fec(order, 3);
  FiniteElementSpace nd_fespace(*mesh, &nd_fec), rt_fespace(*mesh, &rt_fec),
      h1_fespace(*mesh, &h1_fec);

  GridFunction E(nd_fespace, true), B(rt_fespace, true);
  mfem::VectorFunctionCoefficient fer(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::sin(x(1)) + x(2) * x(2);
                                        v(1) = std::cos(x(2)) + x(0);
                                        v(2) = x(0) * x(1) + 1.0;
                                      });
  mfem::VectorFunctionCoefficient fei(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) * x(2) - 0.5;
                                        v(1) = std::sin(x(0)) - x(2);
                                        v(2) = std::cos(x(1)) + x(0) * x(0);
                                      });
  mfem::VectorFunctionCoefficient fbr(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) - 0.3 * x(2);
                                        v(1) = std::sin(x(2)) + 0.5;
                                        v(2) = std::cos(x(0)) - x(1) * x(2);
                                      });
  mfem::VectorFunctionCoefficient fbi(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::cos(x(2)) - 0.2;
                                        v(1) = x(0) * x(2) + 0.1;
                                        v(2) = std::sin(x(1)) - x(0);
                                      });
  E.Real().ProjectCoefficient(fer);
  E.Imag().ProjectCoefficient(fei);
  B.Real().ProjectCoefficient(fbr);
  B.Imag().ProjectCoefficient(fbi);

  const int bdr_attr_max = pmesh.bdr_attributes.Max();
  mfem::Array<int> marker(bdr_attr_max);
  marker = 0;
  marker[7 - 1] = 1;
  const double t_i = 2.0e-3, epsilon_i = 10.0;

  auto TimeIt = [&](auto &&f, int reps)
  {
    // Warm up (also triggers any JIT compilation on GPU backends).
    f();
    double t_min = 1.0e30, t_sum = 0.0;
    for (int r = 0; r < reps; r++)
    {
      auto t0 = std::chrono::steady_clock::now();
      f();
      auto t1 = std::chrono::steady_clock::now();
      const double t = std::chrono::duration<double, std::milli>(t1 - t0).count();
      t_min = std::min(t_min, t);
      t_sum += t;
    }
    return std::make_pair(t_min, t_sum / reps);
  };

  std::cout << "\n=== SurfaceFunctional benchmark: N = " << N << " (" << pmesh.GetNE()
            << " elements, " << N * N << " interface boundary elements, ND order " << order
            << ", " << nd_fespace.GlobalTrueVSize() << " ND dofs) ===\n";

  // 1. Interface dielectric EPR (MA type) over the interior interface.
  {
    volatile double sink;
    auto legacy = [&]()
    {
      // Replicates SurfacePostOperator::GetInterfaceElectricFieldEnergy +
      // GetLocalSurfaceIntegral.
      InterfaceDielectricCoefficient<InterfaceDielectric::MA> f(E, mat_op, t_i, epsilon_i);
      mfem::LinearForm s(&h1_fespace.Get());
      s.AddBoundaryIntegrator(new BoundaryLFIntegrator(f),
                              const_cast<mfem::Array<int> &>(marker));
      s.UseFastAssembly(false);
      s.UseDevice(false);
      s.Assemble();
      s.UseDevice(true);
      double dot = linalg::LocalSum(s);
      Mpi::GlobalSum(1, &dot, comm);
      sink = dot;
    };
    auto t0 = std::chrono::steady_clock::now();
    SurfaceFunctional epr(*mesh, marker, nd_fespace, mat_op, InterfaceDielectric::MA, t_i,
                          epsilon_i);
    auto t1 = std::chrono::steady_clock::now();
    const double t_constr = std::chrono::duration<double, std::milli>(t1 - t0).count();
    auto ceed_path = [&]() { sink = epr.Eval(E); };

    // Verify agreement before timing.
    InterfaceDielectricCoefficient<InterfaceDielectric::MA> f(E, mat_op, t_i, epsilon_i);
    const double ref = RefSurfaceCoefficientIntegral(pmesh, f, marker);
    const double val = epr.Eval(E);
    REQUIRE(val == Catch::Approx(ref).epsilon(1.0e-10));

    auto [legacy_min, legacy_avg] = TimeIt(legacy, n_reps);
    auto [ceed_min, ceed_avg] = TimeIt(ceed_path, n_reps);
    std::cout << "Interface EPR (MA), complex E, per measurement:\n"
              << "  legacy coefficient path: min " << legacy_min << " ms, avg "
              << legacy_avg << " ms\n"
              << "  libCEED functional:      min " << ceed_min << " ms, avg " << ceed_avg
              << " ms  (construction, once: " << t_constr << " ms)\n"
              << "  speedup (avg): " << legacy_avg / ceed_avg << "x\n";
  }

  // 2. Port power over the interior interface (driven solver inner loop quantity).
  {
    volatile double sink;
    mfem::Array<int> attr_list(1);
    attr_list[0] = 7;
    auto legacy = [&]()
    {
      // Replicates LumpedPortData::GetPower.
      std::complex<double> dot;
      {
        RestrictedVectorCoefficient<BdrSurfaceCurrentVectorCoefficient> fbr_c(
            attr_list, B.Real(), mat_op);
        mfem::LinearForm pr(&nd_fespace.Get());
        pr.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fbr_c),
                                 const_cast<mfem::Array<int> &>(marker));
        pr.UseFastAssembly(false);
        pr.UseDevice(false);
        pr.Assemble();
        pr.UseDevice(true);
        dot = std::complex<double>(-(pr * E.Real()), -(pr * E.Imag()));
      }
      {
        RestrictedVectorCoefficient<BdrSurfaceCurrentVectorCoefficient> fbi_c(
            attr_list, B.Imag(), mat_op);
        mfem::LinearForm pi(&nd_fespace.Get());
        pi.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fbi_c),
                                 const_cast<mfem::Array<int> &>(marker));
        pi.UseFastAssembly(false);
        pi.UseDevice(false);
        pi.Assemble();
        pi.UseDevice(true);
        dot += std::complex<double>(-(pi * E.Imag()), pi * E.Real());
      }
      Mpi::GlobalSum(1, &dot, comm);
      sink = dot.real();
    };
    mfem::Vector x0(3);
    x0 = 0.0;
    auto t0 = std::chrono::steady_clock::now();
    SurfaceFunctional power(*mesh, marker, &nd_fespace.Get(), &rt_fespace.Get(), mat_op,
                            SurfaceFlux::POWER, /*two_sided*/ true, x0);
    auto t1 = std::chrono::steady_clock::now();
    const double t_constr = std::chrono::duration<double, std::milli>(t1 - t0).count();
    auto ceed_path = [&]() { sink = power.EvalComplexPower(E, B).real(); };

    auto [legacy_min, legacy_avg] = TimeIt(legacy, n_reps);
    auto [ceed_min, ceed_avg] = TimeIt(ceed_path, n_reps);
    std::cout << "Port power, complex E and B, per measurement:\n"
              << "  legacy linear form path: min " << legacy_min << " ms, avg "
              << legacy_avg << " ms\n"
              << "  libCEED functional:      min " << ceed_min << " ms, avg " << ceed_avg
              << " ms  (construction, once: " << t_constr << " ms)\n"
              << "  speedup (avg): " << legacy_avg / ceed_avg << "x\n";
  }
}

TEST_CASE("SurfaceFunctional nonconformal hex mapped face rules, depths 1 and 2",
          "[surfacefunctional][Serial][GPU]")
{
  // A single order-one H(curl) evaluation keeps the deterministic depth-2 (16-leaf
  // quadrilateral) regression small while retaining a directly comparable depth-1 case.
  MPI_Comm comm = MPI_COMM_WORLD;
  const int refinement_depth = GENERATE(1, 2);
  CAPTURE(refinement_depth);
  fem::DefaultIntegrationOrder::p_trial = 1;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  auto mesh = MakeNCInterfaceMesh(comm, mfem::Element::HEXAHEDRON, refinement_depth);
  auto &pmesh = mesh->Get();
  REQUIRE(pmesh.Nonconforming());

  mfem::ND_FECollection nd_fec(1, 3);
  FiniteElementSpace nd_fespace(*mesh, &nd_fec);
  mfem::ParGridFunction E(&nd_fespace.Get());
  mfem::VectorFunctionCoefficient field(3,
                                        [](const mfem::Vector &x, mfem::Vector &v)
                                        {
                                          v(0) = std::sin(x(1)) + x(2) * x(2);
                                          v(1) = std::cos(x(2)) + x(0);
                                          v(2) = x(0) * x(1) + 1.0;
                                        });
  E.ProjectCoefficient(field);

  mfem::Array<int> marker(pmesh.bdr_attributes.Max());
  marker = 1;
  SurfaceFunctional norm2(SurfaceFunctional::Kind::HCURL_NORM2, *mesh, marker,
                          &nd_fespace.Get());
  const double ref = RefSurfaceHCurlNorm2(pmesh, E, marker);
  CHECK(norm2.Eval(&E) == Catch::Approx(ref).epsilon(1.0e-10));
}

TEST_CASE("SurfaceFunctional Nonconformal Mesh", "[surfacefunctional][Serial][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto elem_type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  auto order = GENERATE(1, 2);
  CAPTURE(elem_type, order);
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  auto mesh = MakeNCInterfaceMesh(comm, elem_type);
  auto &pmesh = mesh->Get();
  REQUIRE(pmesh.Nonconforming());

  config::MaterialData vacuum, dielectric;
  vacuum.attributes = {1};
  dielectric.attributes = {2};
  dielectric.epsilon_r.s[0] = 11.7;
  dielectric.epsilon_r.s[1] = 11.7;
  dielectric.epsilon_r.s[2] = 11.7;
  dielectric.mu_r.s[0] = 1.4;
  dielectric.mu_r.s[1] = 1.4;
  dielectric.mu_r.s[2] = 1.4;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({vacuum, dielectric}, periodic, ProblemType::DRIVEN, *mesh);

  mfem::ND_FECollection nd_fec(order, 3);
  mfem::RT_FECollection rt_fec(order - 1, 3);
  FiniteElementSpace nd_fespace(*mesh, &nd_fec), rt_fespace(*mesh, &rt_fec);

  GridFunction E(nd_fespace, true), B(rt_fespace, true);
  mfem::VectorFunctionCoefficient fer(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::sin(x(1)) + x(2) * x(2);
                                        v(1) = std::cos(x(2)) + x(0);
                                        v(2) = x(0) * x(1) + 1.0;
                                      });
  mfem::VectorFunctionCoefficient fei(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) * x(2) - 0.5;
                                        v(1) = std::sin(x(0)) - x(2);
                                        v(2) = std::cos(x(1)) + x(0) * x(0);
                                      });
  mfem::VectorFunctionCoefficient fbr(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) - 0.3 * x(2);
                                        v(1) = std::sin(x(2)) + 0.5;
                                        v(2) = std::cos(x(0)) - x(1) * x(2);
                                      });
  mfem::VectorFunctionCoefficient fbi(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::cos(x(2)) - 0.2;
                                        v(1) = x(0) * x(2) + 0.1;
                                        v(2) = std::sin(x(1)) - x(0);
                                      });
  E.Real().ProjectCoefficient(fer);
  E.Imag().ProjectCoefficient(fei);
  B.Real().ProjectCoefficient(fbr);
  B.Imag().ProjectCoefficient(fbi);

  const double t_i = 2.0e-3, epsilon_i = 10.0;
  const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
  mfem::Vector x0(3);
  x0(0) = 0.4;
  x0(1) = 0.6;
  x0(2) = -0.2;

  auto CheckPart = [](double v, double r)
  {
    if (std::abs(r) > 1.0e-12)
    {
      CHECK(v == Catch::Approx(r).epsilon(1.0e-10));
    }
    else
    {
      CHECK(std::abs(v) < 1.0e-10);
    }
  };

  // The exterior boundary attribute 1 (z = 0 in the refined half and beyond) and the
  // interior interface attribute 7 both contain a mix of conformal and nonconformal
  // (slave) faces after the partial refinement.
  for (int attr : {1, 6, 7})
  {
    CAPTURE(attr);
    mfem::Array<int> marker(bdr_attr_max);
    marker = 0;
    marker[attr - 1] = 1;

    // Interface dielectric (all variants).
    for (auto type : {InterfaceDielectric::DEFAULT, InterfaceDielectric::MA,
                      InterfaceDielectric::MS, InterfaceDielectric::SA})
    {
      CAPTURE(static_cast<int>(type));
      SurfaceFunctional epr(*mesh, marker, nd_fespace, mat_op, type, t_i, epsilon_i);
      REQUIRE(epr.IsValid());
      auto MakeLegacy = [&]() -> std::unique_ptr<mfem::Coefficient>
      {
        switch (type)
        {
          case InterfaceDielectric::DEFAULT:
            return std::make_unique<
                InterfaceDielectricCoefficient<InterfaceDielectric::DEFAULT>>(
                E, mat_op, t_i, epsilon_i);
          case InterfaceDielectric::MA:
            return std::make_unique<
                InterfaceDielectricCoefficient<InterfaceDielectric::MA>>(E, mat_op, t_i,
                                                                         epsilon_i);
          case InterfaceDielectric::MS:
            return std::make_unique<
                InterfaceDielectricCoefficient<InterfaceDielectric::MS>>(E, mat_op, t_i,
                                                                         epsilon_i);
          case InterfaceDielectric::SA:
            return std::make_unique<
                InterfaceDielectricCoefficient<InterfaceDielectric::SA>>(E, mat_op, t_i,
                                                                         epsilon_i);
        }
        return {};
      };
      auto legacy = MakeLegacy();
      const double ref = RefSurfaceCoefficientIntegral(pmesh, *legacy, marker);
      const double val = epr.Eval(E);
      CAPTURE(ref, val);
      CheckPart(val, ref);
    }

    // Surface flux (all types, both two_sided settings).
    for (auto type : {SurfaceFlux::ELECTRIC, SurfaceFlux::MAGNETIC, SurfaceFlux::POWER})
    {
      for (bool two_sided : {false, true})
      {
        CAPTURE(static_cast<int>(type), two_sided);
        const GridFunction *E_use =
            (type == SurfaceFlux::ELECTRIC || type == SurfaceFlux::POWER) ? &E : nullptr;
        const GridFunction *B_use =
            (type == SurfaceFlux::MAGNETIC || type == SurfaceFlux::POWER) ? &B : nullptr;
        SurfaceFunctional flux(*mesh, marker, E_use ? &nd_fespace.Get() : nullptr,
                               B_use ? &rt_fespace.Get() : nullptr, mat_op, type, two_sided,
                               x0);
        REQUIRE(flux.IsValid());
        auto MakeLegacy =
            [&](const mfem::ParGridFunction *Er,
                const mfem::ParGridFunction *Br) -> std::unique_ptr<mfem::Coefficient>
        {
          switch (type)
          {
            case SurfaceFlux::ELECTRIC:
              return std::make_unique<BdrSurfaceFluxCoefficient<SurfaceFlux::ELECTRIC>>(
                  Er, Br, mat_op, two_sided, x0);
            case SurfaceFlux::MAGNETIC:
              return std::make_unique<BdrSurfaceFluxCoefficient<SurfaceFlux::MAGNETIC>>(
                  Er, Br, mat_op, two_sided, x0);
            case SurfaceFlux::POWER:
              return std::make_unique<BdrSurfaceFluxCoefficient<SurfaceFlux::POWER>>(
                  Er, Br, mat_op, two_sided, x0);
          }
          return {};
        };
        auto legacy_re =
            MakeLegacy(E_use ? &E.Real() : nullptr, B_use ? &B.Real() : nullptr);
        std::complex<double> ref(RefSurfaceCoefficientIntegral(pmesh, *legacy_re, marker),
                                 0.0);
        auto legacy_im =
            MakeLegacy(E_use ? &E.Imag() : nullptr, B_use ? &B.Imag() : nullptr);
        const double ref_im = RefSurfaceCoefficientIntegral(pmesh, *legacy_im, marker);
        if (type == SurfaceFlux::POWER)
        {
          ref += ref_im;
        }
        else
        {
          ref.imag(ref_im);
        }
        const std::complex<double> val = flux.EvalFlux(E_use, B_use);
        CAPTURE(ref.real(), ref.imag(), val.real(), val.imag());
        CheckPart(val.real(), ref.real());
        CheckPart(val.imag(), ref.imag());
      }
    }
  }
}

TEST_CASE("SurfaceFunctional Nonconformal Parallel", "[surfacefunctional][Parallel][GPU]")
{
  // On parallel (possibly nonconformal) meshes, two-sided evaluation across process
  // boundaries uses FaceNbrFieldExchange to pull remote-side physical field values.
  // This test makes that processor-boundary behavior a committed regression: exterior
  // and interior interface surfaces must all assemble through libCEED and match the
  // legacy coefficient oracles.
  MPI_Comm comm = MPI_COMM_WORLD;
  auto elem_type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  const int order = 2;
  CAPTURE(elem_type);
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  auto mesh = MakeNCInterfaceMesh(comm, elem_type);
  auto &pmesh = mesh->Get();

  config::MaterialData vacuum, dielectric;
  vacuum.attributes = {1};
  dielectric.attributes = {2};
  dielectric.epsilon_r.s[0] = 11.7;
  dielectric.epsilon_r.s[1] = 11.7;
  dielectric.epsilon_r.s[2] = 11.7;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({vacuum, dielectric}, periodic, ProblemType::DRIVEN, *mesh);

  mfem::ND_FECollection nd_fec(order, 3);
  mfem::RT_FECollection rt_fec(order - 1, 3);
  FiniteElementSpace nd_fespace(*mesh, &nd_fec), rt_fespace(*mesh, &rt_fec);

  GridFunction E(nd_fespace, true), B(rt_fespace, true);
  mfem::VectorFunctionCoefficient fer(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::sin(x(1)) + x(2) * x(2);
                                        v(1) = std::cos(x(2)) + x(0);
                                        v(2) = x(0) * x(1) + 1.0;
                                      });
  mfem::VectorFunctionCoefficient fei(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) * x(2) - 0.5;
                                        v(1) = std::sin(x(0)) - x(2);
                                        v(2) = std::cos(x(1)) + x(0) * x(0);
                                      });
  mfem::VectorFunctionCoefficient fbr(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) - 0.3 * x(2);
                                        v(1) = std::sin(x(2)) + 0.5;
                                        v(2) = std::cos(x(0)) - x(1) * x(2);
                                      });
  mfem::VectorFunctionCoefficient fbi(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::cos(x(2)) - 0.2;
                                        v(1) = x(0) * x(2) + 0.1;
                                        v(2) = std::sin(x(1)) - x(0);
                                      });
  E.Real().ProjectCoefficient(fer);
  E.Imag().ProjectCoefficient(fei);
  B.Real().ProjectCoefficient(fbr);
  B.Imag().ProjectCoefficient(fbi);

  // The legacy reference evaluates fields on remote sides of shared faces via face
  // neighbor data.
  E.Real().ExchangeFaceNbrData();
  E.Imag().ExchangeFaceNbrData();
  B.Real().ExchangeFaceNbrData();
  B.Imag().ExchangeFaceNbrData();

  const double t_i = 2.0e-3, epsilon_i = 10.0;
  const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;

  for (int attr : {1, 6, 7})
  {
    CAPTURE(attr);
    mfem::Array<int> marker(bdr_attr_max);
    marker = 0;
    marker[attr - 1] = 1;
    auto CheckPart = [](double val, double ref)
    {
      if (std::abs(ref) > 1.0e-12)
      {
        CHECK(val == Catch::Approx(ref).epsilon(1.0e-10));
      }
      else
      {
        CHECK(std::abs(val) < 1.0e-10);
      }
    };

    for (auto type : {InterfaceDielectric::DEFAULT, InterfaceDielectric::MA})
    {
      CAPTURE(static_cast<int>(type));
      SurfaceFunctional epr(*mesh, marker, nd_fespace, mat_op, type, t_i, epsilon_i);

      bool valid = epr.IsValid();
      bool valid_and = valid, valid_or = valid;
      Mpi::GlobalAnd(1, &valid_and, comm);
      Mpi::GlobalOr(1, &valid_or, comm);
      REQUIRE(valid_and == valid_or);
      REQUIRE(valid);

      auto MakeLegacy = [&]() -> std::unique_ptr<mfem::Coefficient>
      {
        if (type == InterfaceDielectric::DEFAULT)
        {
          return std::make_unique<
              InterfaceDielectricCoefficient<InterfaceDielectric::DEFAULT>>(E, mat_op, t_i,
                                                                            epsilon_i);
        }
        return std::make_unique<InterfaceDielectricCoefficient<InterfaceDielectric::MA>>(
            E, mat_op, t_i, epsilon_i);
      };
      auto legacy = MakeLegacy();
      const double ref = RefSurfaceCoefficientIntegral(pmesh, *legacy, marker);
      const double val = epr.Eval(E);
      CAPTURE(ref, val);
      CheckPart(val, ref);
    }

    mfem::Vector x0(3);
    x0 = 0.0;
    SurfaceFunctional power(*mesh, marker, &nd_fespace.Get(), &rt_fespace.Get(), mat_op,
                            SurfaceFlux::POWER, /*two_sided*/ true, x0);
    bool power_valid = power.IsValid();
    bool power_valid_and = power_valid, power_valid_or = power_valid;
    Mpi::GlobalAnd(1, &power_valid_and, comm);
    Mpi::GlobalOr(1, &power_valid_or, comm);
    REQUIRE(power_valid_and == power_valid_or);
    REQUIRE(power_valid);

    auto MakePowerLegacy =
        [&](const mfem::ParGridFunction &Er, const mfem::ParGridFunction &Br)
    {
      auto coeff = std::make_unique<BdrSurfaceFluxCoefficient<SurfaceFlux::POWER>>(
          &Er, &Br, mat_op, /*two_sided*/ true, x0);
      return RefSurfaceCoefficientIntegral(pmesh, *coeff, marker);
    };
    const std::complex<double> ref(
        MakePowerLegacy(E.Real(), B.Real()) + MakePowerLegacy(E.Imag(), B.Imag()), 0.0);
    const std::complex<double> val = power.EvalFlux(&E, &B);
    CAPTURE(ref.real(), val.real(), val.imag());
    CheckPart(val.real(), ref.real());
    CheckPart(val.imag(), 0.0);
  }
}

TEST_CASE("PointFieldEvaluator Domain Fields", "[surfacefunctional][Serial][Parallel][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto elem_type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  auto order = GENERATE(1, 2);
  auto nonconformal = GENERATE(false, true);
  auto complex = GENERATE(false, true);
  CAPTURE(elem_type, order, nonconformal, complex);
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  auto mesh = nonconformal ? MakeNCInterfaceMesh(comm, elem_type)
                           : MakeInterfaceMesh(comm, elem_type);
  auto &pmesh = mesh->Get();

  config::MaterialData vacuum, dielectric;
  vacuum.attributes = {1};
  dielectric.attributes = {2};
  dielectric.epsilon_r.s[0] = 11.7;
  dielectric.epsilon_r.s[1] = 11.7;
  dielectric.epsilon_r.s[2] = 11.7;
  dielectric.mu_r.s[0] = 1.4;
  dielectric.mu_r.s[1] = 1.4;
  dielectric.mu_r.s[2] = 1.4;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({vacuum, dielectric}, periodic, ProblemType::DRIVEN, *mesh);

  mfem::ND_FECollection nd_fec(order, 3);
  mfem::RT_FECollection rt_fec(order - 1, 3);
  FiniteElementSpace nd_fespace(*mesh, &nd_fec), rt_fespace(*mesh, &rt_fec);

  GridFunction E(nd_fespace, complex), B(rt_fespace, complex);
  mfem::VectorFunctionCoefficient fer(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::sin(x(1)) + x(2) * x(2);
                                        v(1) = std::cos(x(2)) + x(0);
                                        v(2) = x(0) * x(1) + 1.0;
                                      });
  mfem::VectorFunctionCoefficient fbr(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) - 0.3 * x(2);
                                        v(1) = std::sin(x(2)) + 0.5;
                                        v(2) = std::cos(x(0)) - x(1) * x(2);
                                      });
  E.Real().ProjectCoefficient(fer);
  B.Real().ProjectCoefficient(fbr);
  if (complex)
  {
    mfem::VectorFunctionCoefficient fei(3,
                                        [](const mfem::Vector &x, mfem::Vector &v)
                                        {
                                          v(0) = x(1) * x(2) - 0.5;
                                          v(1) = std::sin(x(0)) - x(2);
                                          v(2) = std::cos(x(1)) + x(0) * x(0);
                                        });
    mfem::VectorFunctionCoefficient fbi(3,
                                        [](const mfem::Vector &x, mfem::Vector &v)
                                        {
                                          v(0) = std::cos(x(2)) - 0.2;
                                          v(1) = x(0) * x(2) + 0.1;
                                          v(2) = std::sin(x(1)) - x(0);
                                        });
    E.Imag().ProjectCoefficient(fei);
    B.Imag().ProjectCoefficient(fbi);
  }

  // Interpolatory L2 output spaces (legacy ProjectCoefficient on these spaces evaluates
  // the coefficient at the nodal points, which is exactly the libCEED evaluator
  // semantics, at any order).
  mfem::L2_FECollection viz_fec(order, 3);
  mfem::ParFiniteElementSpace viz_scalar(&pmesh, &viz_fec), viz_vector(&pmesh, &viz_fec, 3);

  const double scaling = 2.5;
  auto CheckField = [](const mfem::ParGridFunction &val, const mfem::ParGridFunction &ref)
  {
    // HostRead to sync from device (the evaluator fills on device).
    const double *v = val.HostRead();
    const double *r = ref.HostRead();
    double max_diff = 0.0, max_ref = 0.0;
    for (int i = 0; i < ref.Size(); i++)
    {
      max_diff = std::max(max_diff, std::abs(v[i] - r[i]));
      max_ref = std::max(max_ref, std::abs(r[i]));
    }
    CAPTURE(max_diff, max_ref);
    CHECK(max_diff <= 1.0e-11 * std::max(max_ref, 1.0));
  };

  SECTION("Electric energy density")
  {
    PointFieldEvaluator eval(PointFieldEvaluator::Kind::ENERGY_E, *mesh, mat_op,
                             E.ParFESpace(), nullptr, viz_scalar, scaling);
    REQUIRE(eval.IsValid());
    mfem::ParGridFunction val(&viz_scalar), ref(&viz_scalar);
    eval.Eval(&E, nullptr, val);
    EnergyDensityCoefficient<EnergyDensityType::ELECTRIC> legacy(E, mat_op, scaling);
    ref.ProjectCoefficient(legacy);
    CheckField(val, ref);
  }

  SECTION("Magnetic energy density")
  {
    PointFieldEvaluator eval(PointFieldEvaluator::Kind::ENERGY_M, *mesh, mat_op, nullptr,
                             B.ParFESpace(), viz_scalar, scaling);
    REQUIRE(eval.IsValid());
    mfem::ParGridFunction val(&viz_scalar), ref(&viz_scalar);
    eval.Eval(nullptr, &B, val);
    EnergyDensityCoefficient<EnergyDensityType::MAGNETIC> legacy(B, mat_op, scaling);
    ref.ProjectCoefficient(legacy);
    CheckField(val, ref);
  }

  SECTION("Poynting vector")
  {
    PointFieldEvaluator eval(PointFieldEvaluator::Kind::POYNTING, *mesh, mat_op,
                             E.ParFESpace(), B.ParFESpace(), viz_vector, scaling);
    REQUIRE(eval.IsValid());
    mfem::ParGridFunction val(&viz_vector), ref(&viz_vector);
    eval.Eval(&E, &B, val);
    PoyntingVectorCoefficient legacy(E, B, mat_op, scaling);
    ref.ProjectCoefficient(legacy);
    CheckField(val, ref);
  }
}

TEST_CASE("PointFieldEvaluator Domain Fields 2D",
          "[surfacefunctional][Serial][Parallel][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto elem_type = GENERATE(mfem::Element::TRIANGLE, mfem::Element::QUADRILATERAL);
  auto order = GENERATE(1, 2);
  auto complex = GENERATE(false, true);
  CAPTURE(elem_type, order, complex);
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  auto smesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(3, 2, elem_type, false, 1.2, 0.7));
  smesh->EnsureNodes();
  REQUIRE(Mpi::Size(comm) <= smesh->GetNE());
  auto pmesh_ptr = std::make_unique<mfem::ParMesh>(comm, *smesh);
  Mesh mesh(std::move(pmesh_ptr));
  auto &pmesh = mesh.Get();

  config::MaterialData material;
  material.attributes = {1};
  material.epsilon_r.s[0] = 2.0;
  material.epsilon_r.s[1] = 3.0;
  material.epsilon_r.s[2] = 4.0;
  material.mu_r.s[0] = 1.0;
  material.mu_r.s[1] = 1.5;
  material.mu_r.s[2] = 2.0;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({material}, periodic, ProblemType::DRIVEN, mesh);

  mfem::ND_FECollection nd_fec(order, 2);
  mfem::L2_FECollection l2_fec(order - 1, 2, mfem::BasisType::GaussLegendre,
                               mfem::FiniteElement::INTEGRAL);
  FiniteElementSpace nd_fespace(mesh, &nd_fec), l2_fespace(mesh, &l2_fec);

  GridFunction E(nd_fespace, complex), B(l2_fespace, complex);
  mfem::VectorFunctionCoefficient fer(2,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::sin(x(1)) + 0.2 * x(0);
                                        v(1) = std::cos(x(0)) + x(1) * x(1);
                                      });
  mfem::FunctionCoefficient fbr([](const mfem::Vector &x)
                                { return 0.3 + x(0) - 0.7 * x(1); });
  E.Real().ProjectCoefficient(fer);
  B.Real().ProjectCoefficient(fbr);
  if (complex)
  {
    mfem::VectorFunctionCoefficient fei(2,
                                        [](const mfem::Vector &x, mfem::Vector &v)
                                        {
                                          v(0) = x(0) * x(1) - 0.1;
                                          v(1) = std::sin(x(0) + x(1));
                                        });
    mfem::FunctionCoefficient fbi([](const mfem::Vector &x)
                                  { return std::cos(x(0)) + 0.4 * x(1); });
    E.Imag().ProjectCoefficient(fei);
    B.Imag().ProjectCoefficient(fbi);
  }

  mfem::L2_FECollection viz_fec(order, 2);
  mfem::ParFiniteElementSpace viz_scalar(&pmesh, &viz_fec), viz_vector(&pmesh, &viz_fec, 2);

  const double scaling = 1.7;
  auto CheckField = [](const mfem::ParGridFunction &val, const mfem::ParGridFunction &ref)
  {
    const double *v = val.HostRead();
    const double *r = ref.HostRead();
    double max_diff = 0.0, max_ref = 0.0;
    for (int i = 0; i < ref.Size(); i++)
    {
      max_diff = std::max(max_diff, std::abs(v[i] - r[i]));
      max_ref = std::max(max_ref, std::abs(r[i]));
    }
    CAPTURE(max_diff, max_ref);
    CHECK(max_diff <= 1.0e-11 * std::max(max_ref, 1.0));
  };

  SECTION("Electric energy density")
  {
    PointFieldEvaluator eval(PointFieldEvaluator::Kind::ENERGY_E, mesh, mat_op,
                             E.ParFESpace(), nullptr, viz_scalar, scaling);
    REQUIRE(eval.IsValid());
    mfem::ParGridFunction val(&viz_scalar), ref(&viz_scalar);
    eval.Eval(&E, nullptr, val);
    EnergyDensityCoefficient<EnergyDensityType::ELECTRIC> legacy(E, mat_op, scaling);
    ref.ProjectCoefficient(legacy);
    CheckField(val, ref);
  }

  SECTION("Magnetic flux density")
  {
    PointFieldEvaluator eval(PointFieldEvaluator::Kind::FIELD_B, mesh, mat_op, nullptr,
                             B.ParFESpace(), viz_scalar, 1.0);
    REQUIRE(eval.IsValid());
    mfem::ParGridFunction val(&viz_scalar), ref(&viz_scalar);
    eval.Eval(nullptr, &B, val);
    mfem::GridFunctionCoefficient legacy_real(&B.Real());
    ref.ProjectCoefficient(legacy_real);
    if (complex)
    {
      mfem::ParGridFunction ref_imag(&viz_scalar);
      mfem::GridFunctionCoefficient legacy_imag(&B.Imag());
      ref_imag.ProjectCoefficient(legacy_imag);
      ref += ref_imag;
    }
    CheckField(val, ref);
  }

  SECTION("Magnetic energy density")
  {
    PointFieldEvaluator eval(PointFieldEvaluator::Kind::ENERGY_M, mesh, mat_op, nullptr,
                             B.ParFESpace(), viz_scalar, scaling);
    REQUIRE(eval.IsValid());
    mfem::ParGridFunction val(&viz_scalar), ref(&viz_scalar);
    eval.Eval(nullptr, &B, val);
    EnergyDensityCoefficient<EnergyDensityType::MAGNETIC> legacy(B, mat_op, scaling);
    ref.ProjectCoefficient(legacy);
    CheckField(val, ref);
  }

  SECTION("Poynting vector")
  {
    PointFieldEvaluator eval(PointFieldEvaluator::Kind::POYNTING, mesh, mat_op,
                             E.ParFESpace(), B.ParFESpace(), viz_vector, scaling);
    REQUIRE(eval.IsValid());
    mfem::ParGridFunction val(&viz_vector), ref(&viz_vector);
    eval.Eval(&E, &B, val);
    PoyntingVectorCoefficient legacy(E, B, mat_op, scaling);
    ref.ProjectCoefficient(legacy);
    CheckField(val, ref);
  }
}

TEST_CASE("SurfaceFunctional FarField", "[surfacefunctional][Serial][Parallel][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto elem_type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  auto order = GENERATE(1, 2);
  CAPTURE(elem_type, order);
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  auto mesh = MakeInterfaceMesh(comm, elem_type);
  auto &pmesh = mesh->Get();

  config::MaterialData vacuum, dielectric;
  vacuum.attributes = {1};
  dielectric.attributes = {2};  // Isotropic (far-field requirement)
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({vacuum, dielectric}, periodic, ProblemType::DRIVEN, *mesh);

  mfem::ND_FECollection nd_fec(order, 3);
  mfem::RT_FECollection rt_fec(order - 1, 3);
  FiniteElementSpace nd_fespace(*mesh, &nd_fec), rt_fespace(*mesh, &rt_fec);

  GridFunction E(nd_fespace, true), B(rt_fespace, true);
  mfem::VectorFunctionCoefficient fer(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::sin(x(1)) + x(2) * x(2);
                                        v(1) = std::cos(x(2)) + x(0);
                                        v(2) = x(0) * x(1) + 1.0;
                                      });
  mfem::VectorFunctionCoefficient fei(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) * x(2) - 0.5;
                                        v(1) = std::sin(x(0)) - x(2);
                                        v(2) = std::cos(x(1)) + x(0) * x(0);
                                      });
  mfem::VectorFunctionCoefficient fbr(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) - 0.3 * x(2);
                                        v(1) = std::sin(x(2)) + 0.5;
                                        v(2) = std::cos(x(0)) - x(1) * x(2);
                                      });
  mfem::VectorFunctionCoefficient fbi(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::cos(x(2)) - 0.2;
                                        v(1) = x(0) * x(2) + 0.1;
                                        v(2) = std::sin(x(1)) - x(0);
                                      });
  E.Real().ProjectCoefficient(fer);
  E.Imag().ProjectCoefficient(fei);
  B.Real().ProjectCoefficient(fbr);
  B.Imag().ProjectCoefficient(fbi);

  // Observation directions and (complex) frequency.
  std::vector<std::array<double, 3>> r_naughts;
  for (auto [theta, phi] : {std::pair{0.3, 0.7}, {1.2, 2.1}, {2.4, 4.5}})
  {
    r_naughts.push_back({std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi),
                         std::cos(theta)});
  }
  const double omega_re = 2.7, omega_im = 0.15;

  const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
  mfem::Array<int> marker(bdr_attr_max);
  marker = 0;
  marker[1 - 1] = 1;  // Exterior boundary (z = 0)

  // Legacy reference following GetFarFieldrE.
  std::vector<std::array<double, 3>> integrals_r(r_naughts.size()),
      integrals_i(r_naughts.size());
  for (int i = 0; i < pmesh.GetNBE(); i++)
  {
    if (!marker[pmesh.GetBdrAttribute(i) - 1])
    {
      continue;
    }
    auto *T = const_cast<mfem::ParMesh &>(pmesh).GetBdrElementTransformation(i);
    const auto *fe = nd_fespace.Get().GetBE(i);
    const auto *ir =
        &mfem::IntRules.Get(fe->GetGeomType(), fem::DefaultIntegrationOrder::Get(*T));
    AddStrattonChuIntegrandAtElement(E, B, mat_op, omega_re, omega_im, r_naughts, *T, *ir,
                                     integrals_r, integrals_i);
  }
  Mpi::GlobalSum(3 * r_naughts.size(), integrals_r.data()->data(), comm);
  Mpi::GlobalSum(3 * r_naughts.size(), integrals_i.data()->data(), comm);

  std::vector<std::array<std::complex<double>, 3>> result;
  {
    SurfaceFunctional farfield(*mesh, marker, nd_fespace, rt_fespace, mat_op, r_naughts);
    REQUIRE(farfield.IsValid());
    result = farfield.EvalFarField(E, B, {omega_re, omega_im});
    REQUIRE(result.size() == r_naughts.size());
    for (std::size_t d = 0; d < r_naughts.size(); d++)
    {
      const auto &r = r_naughts[d];
      const auto &Ir = integrals_r[d];
      const auto &Ii = integrals_i[d];
      const std::array<double, 3> cr = {r[1] * Ir[2] - r[2] * Ir[1],
                                        r[2] * Ir[0] - r[0] * Ir[2],
                                        r[0] * Ir[1] - r[1] * Ir[0]};
      const std::array<double, 3> ci = {r[1] * Ii[2] - r[2] * Ii[1],
                                        r[2] * Ii[0] - r[0] * Ii[2],
                                        r[0] * Ii[1] - r[1] * Ii[0]};
      for (int c = 0; c < 3; c++)
      {
        CAPTURE(d, c, cr[c], ci[c], result[d][c].real(), result[d][c].imag());
        CHECK(result[d][c].real() == Catch::Approx(cr[c]).epsilon(1.0e-10).margin(1.0e-14));
        CHECK(result[d][c].imag() == Catch::Approx(ci[c]).epsilon(1.0e-10).margin(1.0e-14));
      }
    }

    // Changing the frequency must reassemble and still agree (different result).
    auto result2 = farfield.EvalFarField(E, B, {2.0 * omega_re, 0.0});
    CHECK(std::abs(result2[0][0] - result[0][0]) > 0.0);
  }

  // Reconstruct against the same live finite element collections. Mapped integration
  // rules retain application lifetime because MFEM caches DofToQuad by rule pointer.
  SurfaceFunctional rebuilt_farfield(*mesh, marker, nd_fespace, rt_fespace, mat_op,
                                     r_naughts);
  REQUIRE(rebuilt_farfield.IsValid());
  const auto rebuilt_result = rebuilt_farfield.EvalFarField(E, B, {omega_re, omega_im});
  REQUIRE(rebuilt_result.size() == result.size());
  for (std::size_t d = 0; d < result.size(); d++)
  {
    for (int c = 0; c < 3; c++)
    {
      CHECK(rebuilt_result[d][c].real() ==
            Catch::Approx(result[d][c].real()).epsilon(1.0e-12).margin(1.0e-14));
      CHECK(rebuilt_result[d][c].imag() ==
            Catch::Approx(result[d][c].imag()).epsilon(1.0e-12).margin(1.0e-14));
    }
  }

  // Interior far-field surfaces are invalid. Construction must make that decision
  // collectively even when only some ranks own elements of the marked interface.
  marker = 0;
  marker[7 - 1] = 1;
  SurfaceFunctional interior_farfield(*mesh, marker, nd_fespace, rt_fespace, mat_op,
                                      r_naughts);
  REQUIRE_FALSE(interior_farfield.IsValid());
}

TEST_CASE("PointFieldEvaluator Boundary Viz Fields",
          "[surfacefunctional][Serial][Parallel][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto elem_type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  auto order = GENERATE(1, 2);
  auto nonconformal = GENERATE(false, true);
  const int refinement_depth = nonconformal ? GENERATE(1, 2) : 1;
  CAPTURE(elem_type, order, nonconformal, refinement_depth);
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  auto mesh = nonconformal ? MakeNCInterfaceMesh(comm, elem_type, refinement_depth)
                           : MakeInterfaceMesh(comm, elem_type);
  auto &pmesh = mesh->Get();

  config::MaterialData vacuum, dielectric;
  vacuum.attributes = {1};
  dielectric.attributes = {2};
  dielectric.epsilon_r.s[0] = 11.7;
  dielectric.epsilon_r.s[1] = 11.7;
  dielectric.epsilon_r.s[2] = 11.7;
  dielectric.mu_r.s[0] = 1.4;
  dielectric.mu_r.s[1] = 1.4;
  dielectric.mu_r.s[2] = 1.4;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({vacuum, dielectric}, periodic, ProblemType::DRIVEN, *mesh);

  mfem::ND_FECollection nd_fec(order, 3);
  mfem::RT_FECollection rt_fec(order - 1, 3);
  FiniteElementSpace nd_fespace(*mesh, &nd_fec), rt_fespace(*mesh, &rt_fec);

  GridFunction E(nd_fespace, false), B(rt_fespace, false);
  mfem::VectorFunctionCoefficient fer(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::sin(x(1)) + x(2) * x(2);
                                        v(1) = std::cos(x(2)) + x(0);
                                        v(2) = x(0) * x(1) + 1.0;
                                      });
  mfem::VectorFunctionCoefficient fbr(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) - 0.3 * x(2);
                                        v(1) = std::sin(x(2)) + 0.5;
                                        v(2) = std::cos(x(0)) - x(1) * x(2);
                                      });
  E.Real().ProjectCoefficient(fer);
  B.Real().ProjectCoefficient(fbr);
  E.Real().ExchangeFaceNbrData();
  B.Real().ExchangeFaceNbrData();

  const int lod = order;
  const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
  mfem::Array<int> marker(bdr_attr_max);
  marker = 1;
  const long long plan_builds_before = FaceSamplingPlan::BuildCount();
  const long long plan_reuses_before = FaceSamplingPlan::ReuseCount();
  auto sampling_plan = std::make_shared<FaceSamplingPlan>(*mesh, marker, lod);
  REQUIRE(FaceSamplingPlan::BuildCount() == plan_builds_before + 1);

  auto TestKind = [&](PointFieldEvaluator::Kind kind,
                      const mfem::ParFiniteElementSpace &fes,
                      const mfem::ParGridFunction &U)
  {
    PointFieldEvaluator viz(kind, *mesh, marker, fes, lod, sampling_plan);
    REQUIRE(viz.SamplingPlan() == sampling_plan.get());

    // The validity decision must be identical on all ranks; interior surfaces split
    // across processes fall back (consistently).
    bool valid = viz.IsValid();
    bool valid_and = valid, valid_or = valid;
    Mpi::GlobalAnd(1, &valid_and, comm);
    Mpi::GlobalOr(1, &valid_or, comm);
    REQUIRE(valid_and == valid_or);
    REQUIRE(valid);

    Vector buffer(viz.BufferSize());
    buffer.UseDevice(true);
    // Constructor warm-up releases its source storage. Two evaluations exercise both
    // the first direct reattachment and the later take-and-repoint path, including MPI
    // face-neighbor exporters on partition interfaces.
    viz.EvalBuffer(U, buffer);
    viz.EvalBuffer(U, buffer);
    const double *buf = buffer.HostRead();
    const auto &bases = viz.BufferBases();
    const int component_stride = viz.BufferSize() / 3;

    BdrFieldVectorCoefficient legacy(U);
    mfem::Vector V(3);
    for (int i = 0; i < pmesh.GetNBE(); i++)
    {
      const auto &RefG =
          *mfem::GlobGeometryRefiner.Refine(pmesh.GetBdrElementGeometry(i), lod, 1);
      auto *T = pmesh.GetBdrElementTransformation(i);
      for (int j = 0; j < RefG.RefPts.GetNPoints(); j++)
      {
        const auto &ip = RefG.RefPts.IntPoint(j);
        T->SetIntPoint(&ip);
        legacy.Eval(V, *T, ip);
        for (int c = 0; c < 3; c++)
        {
          const double val = buf[bases[i] + j + c * component_stride];
          CAPTURE(i, j, c, V(c), val);
          CHECK(val == Catch::Approx(V(c)).epsilon(1.0e-10).margin(1.0e-13));
        }
      }
    }
  };

  TestKind(PointFieldEvaluator::Kind::FIELD_E, nd_fespace, E.Real());
  TestKind(PointFieldEvaluator::Kind::FIELD_B, rt_fespace, B.Real());

  // Material-dependent boundary visualization kinds (surface charge, surface current,
  // boundary energy densities) against the corresponding legacy coefficients.
  const double scaling = 1.7;
  auto TestKindCoeff = [&](PointFieldEvaluator::Kind kind,
                           const mfem::ParFiniteElementSpace &fes,
                           mfem::Coefficient *legacy_s, mfem::VectorCoefficient *legacy_v)
  {
    PointFieldEvaluator viz(kind, *mesh, marker, fes, mat_op, lod, scaling, sampling_plan);
    REQUIRE(viz.SamplingPlan() == sampling_plan.get());
    const int nc = viz.BufferNumComp();
    bool valid = viz.IsValid();
    bool valid_and = valid, valid_or = valid;
    Mpi::GlobalAnd(1, &valid_and, comm);
    Mpi::GlobalOr(1, &valid_or, comm);
    REQUIRE(valid_and == valid_or);
    REQUIRE(valid);
    Vector buffer(viz.BufferSize());
    buffer.UseDevice(true);
    const auto &u = kind == PointFieldEvaluator::Kind::CURRENT_J ||
                            kind == PointFieldEvaluator::Kind::ENERGY_M
                        ? B.Real()
                        : E.Real();
    viz.EvalBuffer(u, buffer);
    viz.EvalBuffer(u, buffer);
    const double *buf = buffer.HostRead();
    const auto &bases = viz.BufferBases();
    const int component_stride = viz.BufferSize() / nc;
    mfem::Vector V(3);
    for (int i = 0; i < pmesh.GetNBE(); i++)
    {
      const auto &RefG =
          *mfem::GlobGeometryRefiner.Refine(pmesh.GetBdrElementGeometry(i), lod, 1);
      auto *T = pmesh.GetBdrElementTransformation(i);
      for (int j = 0; j < RefG.RefPts.GetNPoints(); j++)
      {
        const auto &ip = RefG.RefPts.IntPoint(j);
        T->SetIntPoint(&ip);
        if (nc == 1)
        {
          const double ref = legacy_s->Eval(*T, ip);
          const double val = buf[bases[i] + j];
          CAPTURE(i, j, ref, val);
          CHECK(val == Catch::Approx(ref).epsilon(1.0e-10).margin(1.0e-13));
        }
        else
        {
          legacy_v->Eval(V, *T, ip);
          for (int c = 0; c < 3; c++)
          {
            const double val = buf[bases[i] + j + c * component_stride];
            CAPTURE(i, j, c, V(c), val);
            CHECK(val == Catch::Approx(V(c)).epsilon(1.0e-10).margin(1.0e-13));
          }
        }
      }
    }
  };
  {
    BdrSurfaceFluxCoefficient<SurfaceFlux::ELECTRIC> q_legacy(
        &E.Real(), nullptr, mat_op, true, mfem::Vector(), scaling);
    TestKindCoeff(PointFieldEvaluator::Kind::FLUX_Q, nd_fespace, &q_legacy, nullptr);
  }
  {
    BdrSurfaceCurrentVectorCoefficient j_legacy(B.Real(), mat_op, scaling);
    TestKindCoeff(PointFieldEvaluator::Kind::CURRENT_J, rt_fespace, nullptr, &j_legacy);
  }
  {
    EnergyDensityCoefficient<EnergyDensityType::ELECTRIC> ue_legacy(E, mat_op, scaling);
    TestKindCoeff(PointFieldEvaluator::Kind::ENERGY_E, nd_fespace, &ue_legacy, nullptr);
  }
  {
    EnergyDensityCoefficient<EnergyDensityType::MAGNETIC> um_legacy(B, mat_op, scaling);
    TestKindCoeff(PointFieldEvaluator::Kind::ENERGY_M, rt_fespace, &um_legacy, nullptr);
  }
  {
    PointFieldEvaluator poynting(PointFieldEvaluator::Kind::POYNTING, *mesh, marker,
                                 nd_fespace, rt_fespace, mat_op, lod, scaling,
                                 sampling_plan);
    REQUIRE(poynting.SamplingPlan() == sampling_plan.get());
    REQUIRE(poynting.IsValid());
    Vector buffer(poynting.BufferSize());
    buffer.UseDevice(true);
    poynting.EvalBuffer(&E, &B, buffer);
  }

  // Trace-backed evaluators consume canonical, side-separated physical ND/RT traces
  // through EVAL_NONE. Compare every independent derived boundary output to the legacy
  // evaluator. The interface mesh exercises two local sides in serial and face-neighbor
  // (ghost) routes in parallel; this test is enabled for the device backends too.
  auto trace_cache =
      std::make_shared<BoundaryPhysicalTraceCache>(*mesh, marker, sampling_plan);
  auto CheckTrace =
      [&](PointFieldEvaluator &legacy, PointFieldEvaluator &traced, const Vector &u)
  {
    Vector expected(legacy.BufferSize()), actual(traced.BufferSize());
    expected.UseDevice(true);
    actual.UseDevice(true);
    legacy.EvalBuffer(u, expected);
    traced.EvalBuffer(u, actual);
    REQUIRE(expected.Size() == actual.Size());
    const double *expected_values = expected.HostRead();
    const double *actual_values = actual.HostRead();
    for (int i = 0; i < expected.Size(); i++)
    {
      CAPTURE(i, expected_values[i], actual_values[i]);
      CHECK(actual_values[i] ==
            Catch::Approx(expected_values[i]).epsilon(1.0e-10).margin(1.0e-13));
    }
  };
  PointFieldEvaluator e_legacy(PointFieldEvaluator::Kind::FIELD_E, *mesh, marker,
                               nd_fespace, lod, sampling_plan);
  PointFieldEvaluator e_traced(PointFieldEvaluator::Kind::FIELD_E, *mesh, marker,
                               nd_fespace, lod, sampling_plan, trace_cache);
  PointFieldEvaluator b_legacy(PointFieldEvaluator::Kind::FIELD_B, *mesh, marker,
                               rt_fespace, lod, sampling_plan);
  PointFieldEvaluator b_traced(PointFieldEvaluator::Kind::FIELD_B, *mesh, marker,
                               rt_fespace, lod, sampling_plan, trace_cache);
  PointFieldEvaluator q_legacy_eval(PointFieldEvaluator::Kind::FLUX_Q, *mesh, marker,
                                    nd_fespace, mat_op, lod, scaling, sampling_plan);
  PointFieldEvaluator q_traced_eval(PointFieldEvaluator::Kind::FLUX_Q, *mesh, marker,
                                    nd_fespace, mat_op, lod, scaling, sampling_plan,
                                    trace_cache);
  PointFieldEvaluator j_legacy_eval(PointFieldEvaluator::Kind::CURRENT_J, *mesh, marker,
                                    rt_fespace, mat_op, lod, scaling, sampling_plan);
  PointFieldEvaluator j_traced_eval(PointFieldEvaluator::Kind::CURRENT_J, *mesh, marker,
                                    rt_fespace, mat_op, lod, scaling, sampling_plan,
                                    trace_cache);
  PointFieldEvaluator ue_legacy_eval(PointFieldEvaluator::Kind::ENERGY_E, *mesh, marker,
                                     nd_fespace, mat_op, lod, scaling, sampling_plan);
  PointFieldEvaluator ue_traced_eval(PointFieldEvaluator::Kind::ENERGY_E, *mesh, marker,
                                     nd_fespace, mat_op, lod, scaling, sampling_plan,
                                     trace_cache);
  PointFieldEvaluator um_legacy_eval(PointFieldEvaluator::Kind::ENERGY_M, *mesh, marker,
                                     rt_fespace, mat_op, lod, scaling, sampling_plan);
  PointFieldEvaluator um_traced_eval(PointFieldEvaluator::Kind::ENERGY_M, *mesh, marker,
                                     rt_fespace, mat_op, lod, scaling, sampling_plan,
                                     trace_cache);
  PointFieldEvaluator s_legacy_eval(PointFieldEvaluator::Kind::POYNTING, *mesh, marker,
                                    nd_fespace, rt_fespace, mat_op, lod, scaling,
                                    sampling_plan);
  PointFieldEvaluator s_traced_eval(PointFieldEvaluator::Kind::POYNTING, *mesh, marker,
                                    nd_fespace, rt_fespace, mat_op, lod, scaling,
                                    sampling_plan, trace_cache);
  REQUIRE(e_traced.TraceCache() == trace_cache.get());
  REQUIRE(b_traced.TraceCache() == trace_cache.get());

  // The fused bundle consumes the two physical trace families once, then exposes only
  // the requested packed slice to each PointFieldEvaluator. Compare every slice against
  // the current independent trace-backed evaluators on exterior/interior and MPI ghost
  // routes exercised by this fixture.
  const long long bundle_count_before = BoundaryDerivedFieldBundle::BundleCount();
  const long long bundle_applies_before = BoundaryDerivedFieldBundle::QFunctionApplyCount();
  auto derived_bundle = std::make_shared<BoundaryDerivedFieldBundle>(
      *mesh, marker, mat_op, nd_fespace, rt_fespace, sampling_plan, trace_cache, E, B,
      scaling, scaling);
  REQUIRE(derived_bundle->IsValid());
  REQUIRE(BoundaryDerivedFieldBundle::BundleCount() == bundle_count_before + 1);
  PointFieldEvaluator e_bundle(PointFieldEvaluator::Kind::FIELD_E, *mesh, marker,
                               nd_fespace, lod, sampling_plan, trace_cache, derived_bundle);
  PointFieldEvaluator b_bundle(PointFieldEvaluator::Kind::FIELD_B, *mesh, marker,
                               rt_fespace, lod, sampling_plan, trace_cache, derived_bundle);
  PointFieldEvaluator q_bundle(PointFieldEvaluator::Kind::FLUX_Q, *mesh, marker, nd_fespace,
                               mat_op, lod, scaling, sampling_plan, trace_cache,
                               derived_bundle);
  PointFieldEvaluator j_bundle(PointFieldEvaluator::Kind::CURRENT_J, *mesh, marker,
                               rt_fespace, mat_op, lod, scaling, sampling_plan, trace_cache,
                               derived_bundle);
  PointFieldEvaluator ue_bundle(PointFieldEvaluator::Kind::ENERGY_E, *mesh, marker,
                                nd_fespace, mat_op, lod, scaling, sampling_plan,
                                trace_cache, derived_bundle);
  PointFieldEvaluator um_bundle(PointFieldEvaluator::Kind::ENERGY_M, *mesh, marker,
                                rt_fespace, mat_op, lod, scaling, sampling_plan,
                                trace_cache, derived_bundle);
  PointFieldEvaluator s_bundle(PointFieldEvaluator::Kind::POYNTING, *mesh, marker,
                               nd_fespace, rt_fespace, mat_op, lod, scaling, sampling_plan,
                               trace_cache, derived_bundle);
  REQUIRE(e_bundle.SamplingPlan() == sampling_plan.get());
  REQUIRE(e_bundle.TraceCache() == trace_cache.get());
  CheckTrace(e_traced, e_bundle, E.Real());
  CheckTrace(b_traced, b_bundle, B.Real());
  CheckTrace(q_traced_eval, q_bundle, E.Real());
  CheckTrace(j_traced_eval, j_bundle, B.Real());
  CheckTrace(ue_traced_eval, ue_bundle, E.Real());
  CheckTrace(um_traced_eval, um_bundle, B.Real());
  {
    Vector expected(s_traced_eval.BufferSize()), actual(s_bundle.BufferSize());
    expected.UseDevice(true);
    actual.UseDevice(true);
    s_traced_eval.EvalBuffer(&E, &B, expected);
    s_bundle.EvalBuffer(&E, &B, actual);
    const double *expected_values = expected.HostRead();
    const double *actual_values = actual.HostRead();
    REQUIRE(expected.Size() == actual.Size());
    for (int i = 0; i < expected.Size(); i++)
    {
      CAPTURE(i, expected_values[i], actual_values[i]);
      CHECK(actual_values[i] ==
            Catch::Approx(expected_values[i]).epsilon(1.0e-10).margin(1.0e-13));
    }
  }
  REQUIRE(BoundaryDerivedFieldBundle::QFunctionApplyCount() > bundle_applies_before);
  REQUIRE(BoundaryDerivedFieldBundle::AvoidedIndependentApplyCount() > 0);
  REQUIRE(BoundaryDerivedFieldBundle::PackedByteCount() > 0);
  const long long bundle_applies_cached = BoundaryDerivedFieldBundle::QFunctionApplyCount();
  Vector cached(e_bundle.BufferSize());
  cached.UseDevice(true);
  e_bundle.EvalBuffer(E.Real(), cached);
  REQUIRE(BoundaryDerivedFieldBundle::QFunctionApplyCount() == bundle_applies_cached);

  // The complex bundle consumes the separately materialized real/imaginary physical
  // traces in one derived QFunction apply per route. Its public linear slices must stay
  // phase-resolved, while U_e/U_m/S retain the old same-phase accumulation exactly.
  GridFunction E_complex(nd_fespace, true), B_complex(rt_fespace, true);
  E_complex.Real().ProjectCoefficient(fer);
  B_complex.Real().ProjectCoefficient(fbr);
  mfem::VectorFunctionCoefficient fei(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(0) - 0.4 * x(1) * x(2);
                                        v(1) = std::sin(x(0) + x(2));
                                        v(2) = 0.25 + x(1) * x(1);
                                      });
  mfem::VectorFunctionCoefficient fbi(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::cos(x(1)) - x(2);
                                        v(1) = x(0) * x(2) + 0.2;
                                        v(2) = std::sin(x(0) - x(1));
                                      });
  E_complex.Imag().ProjectCoefficient(fei);
  B_complex.Imag().ProjectCoefficient(fbi);
  E_complex.Real().ExchangeFaceNbrData();
  E_complex.Imag().ExchangeFaceNbrData();
  B_complex.Real().ExchangeFaceNbrData();
  B_complex.Imag().ExchangeFaceNbrData();

  auto complex_trace_cache =
      std::make_shared<BoundaryPhysicalTraceCache>(*mesh, marker, sampling_plan);
  const long long complex_qfunction_before =
      BoundaryDerivedFieldBundle::QFunctionApplyCount();
  const long long complex_phases_before = BoundaryDerivedFieldBundle::PhaseCount();
  const long long complex_phase_avoided_before =
      BoundaryDerivedFieldBundle::AvoidedPhaseOperatorApplyCount();
  auto complex_bundle = std::make_shared<BoundaryDerivedFieldBundle>(
      *mesh, marker, mat_op, nd_fespace, rt_fespace, sampling_plan, complex_trace_cache,
      E_complex, B_complex, scaling, scaling);
  REQUIRE(complex_bundle->IsValid());
  REQUIRE(complex_bundle->BatchesComplexPhases(E_complex.Real(), B_complex.Real()));
  REQUIRE(complex_bundle->BatchesComplexPhases(E_complex.Imag(), B_complex.Imag()));
  REQUIRE_THROWS(
      complex_bundle->RegisteredSlice(PointFieldKind::ENERGY_E, /*imaginary_phase*/ true));

  PointFieldEvaluator e_complex_traced(PointFieldEvaluator::Kind::FIELD_E, *mesh, marker,
                                       nd_fespace, lod, sampling_plan, complex_trace_cache);
  PointFieldEvaluator b_complex_traced(PointFieldEvaluator::Kind::FIELD_B, *mesh, marker,
                                       rt_fespace, lod, sampling_plan, complex_trace_cache);
  PointFieldEvaluator q_complex_traced(PointFieldEvaluator::Kind::FLUX_Q, *mesh, marker,
                                       nd_fespace, mat_op, lod, scaling, sampling_plan,
                                       complex_trace_cache);
  PointFieldEvaluator j_complex_traced(PointFieldEvaluator::Kind::CURRENT_J, *mesh, marker,
                                       rt_fespace, mat_op, lod, scaling, sampling_plan,
                                       complex_trace_cache);
  PointFieldEvaluator ue_complex_traced(PointFieldEvaluator::Kind::ENERGY_E, *mesh, marker,
                                        nd_fespace, mat_op, lod, scaling, sampling_plan,
                                        complex_trace_cache);
  PointFieldEvaluator um_complex_traced(PointFieldEvaluator::Kind::ENERGY_M, *mesh, marker,
                                        rt_fespace, mat_op, lod, scaling, sampling_plan,
                                        complex_trace_cache);
  PointFieldEvaluator s_complex_traced(PointFieldEvaluator::Kind::POYNTING, *mesh, marker,
                                       nd_fespace, rt_fespace, mat_op, lod, scaling,
                                       sampling_plan, complex_trace_cache);
  PointFieldEvaluator e_complex_bundle(PointFieldEvaluator::Kind::FIELD_E, *mesh, marker,
                                       nd_fespace, lod, sampling_plan, complex_trace_cache,
                                       complex_bundle);
  PointFieldEvaluator b_complex_bundle(PointFieldEvaluator::Kind::FIELD_B, *mesh, marker,
                                       rt_fespace, lod, sampling_plan, complex_trace_cache,
                                       complex_bundle);
  PointFieldEvaluator q_complex_bundle(PointFieldEvaluator::Kind::FLUX_Q, *mesh, marker,
                                       nd_fespace, mat_op, lod, scaling, sampling_plan,
                                       complex_trace_cache, complex_bundle);
  PointFieldEvaluator j_complex_bundle(PointFieldEvaluator::Kind::CURRENT_J, *mesh, marker,
                                       rt_fespace, mat_op, lod, scaling, sampling_plan,
                                       complex_trace_cache, complex_bundle);
  PointFieldEvaluator ue_complex_bundle(PointFieldEvaluator::Kind::ENERGY_E, *mesh, marker,
                                        nd_fespace, mat_op, lod, scaling, sampling_plan,
                                        complex_trace_cache, complex_bundle);
  PointFieldEvaluator um_complex_bundle(PointFieldEvaluator::Kind::ENERGY_M, *mesh, marker,
                                        rt_fespace, mat_op, lod, scaling, sampling_plan,
                                        complex_trace_cache, complex_bundle);
  PointFieldEvaluator s_complex_bundle(PointFieldEvaluator::Kind::POYNTING, *mesh, marker,
                                       nd_fespace, rt_fespace, mat_op, lod, scaling,
                                       sampling_plan, complex_trace_cache, complex_bundle);

  // Each vector source selects the matching linear packed phase slice, not an
  // accumulated real-plus-imaginary field.
  CheckTrace(e_complex_traced, e_complex_bundle, E_complex.Real());
  const long long complex_qfunction_after_real =
      BoundaryDerivedFieldBundle::QFunctionApplyCount();
  CheckTrace(e_complex_traced, e_complex_bundle, E_complex.Imag());
  CheckTrace(b_complex_traced, b_complex_bundle, B_complex.Real());
  CheckTrace(b_complex_traced, b_complex_bundle, B_complex.Imag());
  CheckTrace(q_complex_traced, q_complex_bundle, E_complex.Real());
  CheckTrace(q_complex_traced, q_complex_bundle, E_complex.Imag());
  CheckTrace(j_complex_traced, j_complex_bundle, B_complex.Real());
  CheckTrace(j_complex_traced, j_complex_bundle, B_complex.Imag());
  REQUIRE(BoundaryDerivedFieldBundle::QFunctionApplyCount() ==
          complex_qfunction_after_real);

  auto CheckComplexCombined = [&](PointFieldEvaluator &traced, PointFieldEvaluator &bundled,
                                  const GridFunction *single_source)
  {
    Vector expected(traced.BufferSize()), actual(bundled.BufferSize());
    expected.UseDevice(true);
    actual.UseDevice(true);
    if (single_source)
    {
      traced.EvalBuffer(*single_source, expected);
      bundled.EvalBuffer(*single_source, actual);
    }
    else
    {
      traced.EvalBuffer(&E_complex, &B_complex, expected);
      bundled.EvalBuffer(&E_complex, &B_complex, actual);
    }
    REQUIRE(expected.Size() == actual.Size());
    const double *expected_values = expected.HostRead();
    const double *actual_values = actual.HostRead();
    for (int i = 0; i < expected.Size(); i++)
    {
      CAPTURE(i, expected_values[i], actual_values[i]);
      CHECK(actual_values[i] ==
            Catch::Approx(expected_values[i]).epsilon(1.0e-10).margin(1.0e-13));
    }
  };
  CheckComplexCombined(ue_complex_traced, ue_complex_bundle, &E_complex);
  CheckComplexCombined(um_complex_traced, um_complex_bundle, &B_complex);
  CheckComplexCombined(s_complex_traced, s_complex_bundle, nullptr);

  const long long complex_qfunction_applies =
      BoundaryDerivedFieldBundle::QFunctionApplyCount() - complex_qfunction_before;
  const long long real_route_applies = bundle_applies_cached - bundle_applies_before;
  REQUIRE(complex_qfunction_applies > 0);
  // One complex launch per route replaces the two independent real/imaginary launches
  // used by commit 3's one-phase packing.
  REQUIRE(complex_qfunction_applies == real_route_applies);
  REQUIRE(BoundaryDerivedFieldBundle::PhaseCount() - complex_phases_before ==
          2 * complex_qfunction_applies);
  REQUIRE(BoundaryDerivedFieldBundle::AvoidedPhaseOperatorApplyCount() -
              complex_phase_avoided_before ==
          complex_qfunction_applies);

  CheckTrace(e_legacy, e_traced, E.Real());
  CheckTrace(b_legacy, b_traced, B.Real());
  CheckTrace(q_legacy_eval, q_traced_eval, E.Real());
  CheckTrace(j_legacy_eval, j_traced_eval, B.Real());
  CheckTrace(ue_legacy_eval, ue_traced_eval, E.Real());
  CheckTrace(um_legacy_eval, um_traced_eval, B.Real());
  {
    Vector expected(s_legacy_eval.BufferSize()), actual(s_traced_eval.BufferSize());
    expected.UseDevice(true);
    actual.UseDevice(true);
    s_legacy_eval.EvalBuffer(&E, &B, expected);
    s_traced_eval.EvalBuffer(&E, &B, actual);
    const double *expected_values = expected.HostRead();
    const double *actual_values = actual.HostRead();
    for (int i = 0; i < expected.Size(); i++)
    {
      CAPTURE(i, expected_values[i], actual_values[i]);
      CHECK(actual_values[i] ==
            Catch::Approx(expected_values[i]).epsilon(1.0e-10).margin(1.0e-13));
    }
  }
  REQUIRE(BoundaryPhysicalTraceCache::OperatorApplyCount() >= 4);
  REQUIRE(BoundaryPhysicalTraceCache::CacheHitCount() >= 1);
  // The cache key includes an explicit save generation, not only the Vector address.
  // This models a solver reusing the same storage for its next ParaView save.
  const long long applies_before_invalidate =
      BoundaryPhysicalTraceCache::OperatorApplyCount();
  trace_cache->Invalidate();
  Vector regenerated(e_traced.BufferSize());
  regenerated.UseDevice(true);
  e_traced.EvalBuffer(E.Real(), regenerated);
  REQUIRE(BoundaryPhysicalTraceCache::OperatorApplyCount() > applies_before_invalidate);

  // A direct boundary evaluator retains the existing API and constructs a private plan.
  // Its field output must remain equivalent to an evaluator using the shared plan.
  PointFieldEvaluator private_eval(PointFieldEvaluator::Kind::FIELD_E, *mesh, marker,
                                   nd_fespace, lod);
  PointFieldEvaluator shared_eval(PointFieldEvaluator::Kind::FIELD_E, *mesh, marker,
                                  nd_fespace, lod, sampling_plan);
  REQUIRE(private_eval.SamplingPlan());
  REQUIRE(private_eval.SamplingPlan() != sampling_plan.get());
  REQUIRE(shared_eval.SamplingPlan() == sampling_plan.get());
  REQUIRE(FaceSamplingPlan::BuildCount() == plan_builds_before + 2);
  REQUIRE(FaceSamplingPlan::ReuseCount() >= plan_reuses_before + 7);
  REQUIRE(private_eval.IsValid());
  REQUIRE(shared_eval.IsValid());
  Vector private_buffer(private_eval.BufferSize()), shared_buffer(shared_eval.BufferSize());
  private_buffer.UseDevice(true);
  shared_buffer.UseDevice(true);
  private_eval.EvalBuffer(E.Real(), private_buffer);
  shared_eval.EvalBuffer(E.Real(), shared_buffer);
  REQUIRE(private_buffer.Size() == shared_buffer.Size());
  const double *private_values = private_buffer.HostRead();
  const double *shared_values = shared_buffer.HostRead();
  for (int i = 0; i < private_buffer.Size(); i++)
  {
    CHECK(private_values[i] == Catch::Approx(shared_values[i]).margin(1.0e-13));
  }

  // A plan cannot outlive in-place mesh topology/geometry revisions. NodesUpdated is a
  // cheap deterministic invalidation check; AMR/rebalance also changes GetSequence().
  REQUIRE(sampling_plan->Matches(*mesh, marker, lod));
  pmesh.NodesUpdated();
  REQUIRE_FALSE(sampling_plan->Matches(*mesh, marker, lod));
  // MFEM's fatal-error path is locally catchable in serial but terminates a non-root
  // rank before Catch2 can report it. The stale-plan predicate above is collective
  // evidence in MPI; retain the throwing trace-cache contract check in serial.
  if (Mpi::Size(comm) == 1)
  {
    REQUIRE_THROWS(trace_cache->Get(nd_fespace, E.Real()));
  }
}

TEST_CASE("Boundary trace NC coarse-face unions 2D",
          "[surfacefunctional][Serial][Parallel][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const auto elem_type = GENERATE(mfem::Element::TRIANGLE, mfem::Element::QUADRILATERAL);
  const int refinement_depth = GENERATE(1, 2);
  CAPTURE(elem_type, refinement_depth);
  fem::DefaultIntegrationOrder::p_trial = 2;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  const int coarse_size = Mpi::Size(comm) >= 8 ? 4 : 2;
  CheckBoundaryTraceNCUnion(
      MakeNCInterfaceMesh2D(comm, elem_type, refinement_depth, coarse_size));
}

TEST_CASE("Boundary trace NC coarse-face unions 3D",
          "[surfacefunctional][Serial][Parallel][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const auto elem_type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  const int refinement_depth = GENERATE(1, 2);
  CAPTURE(elem_type, refinement_depth);
  fem::DefaultIntegrationOrder::p_trial = 2;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  const int coarse_size = Mpi::Size(comm) >= 8 ? 3 : 2;
  CheckBoundaryTraceNCUnion(
      MakeNCCoarseInterfaceMesh3D(comm, elem_type, refinement_depth, coarse_size));
}

TEST_CASE("Boundary trace NC mixed wedge/pyramid coarse-owner unions",
          "[surfacefunctional][Serial][Parallel][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const auto interface =
      GENERATE(MixedCoarseInterface::PRISM_QUAD, MixedCoarseInterface::PRISM_TRIANGLE,
               MixedCoarseInterface::PYRAMID_QUAD, MixedCoarseInterface::PYRAMID_TRIANGLE);
  const int refinement_depth = GENERATE(1, 2);
  CAPTURE(static_cast<int>(interface), refinement_depth);
  fem::DefaultIntegrationOrder::p_trial = 2;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  CheckBoundaryTraceNCUnion(MakeNCMixedCoarseOwnerMesh(comm, interface, refinement_depth),
                            CoarseGeometry(interface), CoarseFaceGeometry(interface),
                            Mpi::Size(comm) > 1);
}

TEST_CASE("FaceNbrFieldExchange 2D", "[surfacefunctional][Serial][Parallel][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto mesh = MakeInterfaceMesh2D(comm, mfem::Element::TRIANGLE);
  auto &pmesh = mesh->Get();

  constexpr int order = 2;
  mfem::ND_FECollection nd_fec(order, pmesh.Dimension());
  FiniteElementSpace nd_fespace(*mesh, &nd_fec);
  mfem::ParGridFunction E(&nd_fespace.Get());
  mfem::VectorFunctionCoefficient field(2,
                                        [](const mfem::Vector &x, mfem::Vector &v)
                                        {
                                          v(0) = std::sin(x(1)) + 0.25 * x(0);
                                          v(1) = std::cos(x(0)) + x(0) * x(1);
                                        });
  E.ProjectCoefficient(field);

  const int num_ghost = pmesh.GetNFaceNeighborElements();
  std::vector<FaceNbrFieldExchange::Request> requests;
  for (int fn = 0; fn < num_ghost; fn++)
  {
    auto &req = requests.emplace_back();
    req.face_nbr_elem = fn;
    req.source_mask = 0b01u;
    req.point_key = {mfem::Element::TRIANGLE, order, 3};
    req.pts.resize(3);
    req.pts[0].Set2(0.1, 0.2);
    req.pts[1].Set2(0.2, 0.1);
    req.pts[2].Set2(0.25, 0.25);
  }
  FaceNbrFieldExchange exchange(*mesh, {&nd_fespace.Get(), nullptr, nullptr, nullptr},
                                requests);
  exchange.Exchange({&E, nullptr, nullptr, nullptr});

  E.ExchangeFaceNbrData();
  const double *values = exchange.Imported().HostRead();
  mfem::Vector ref(2);
  int num_checked = 0;
  for (std::size_t r = 0; r < requests.size(); r++)
  {
    const auto &req = requests[r];
    const int offset = exchange.ImportOffset(static_cast<int>(r), 0);
    REQUIRE(offset >= 0);
    for (std::size_t j = 0; j < req.pts.size(); j++)
    {
      E.GetVectorValue(pmesh.GetNE() + req.face_nbr_elem, req.pts[j], ref);
      for (int c = 0; c < 2; c++)
      {
        CAPTURE(r, j, c);
        CHECK(values[offset + 2 * j + c] == Catch::Approx(ref(c)).margin(1.0e-12));
        num_checked++;
      }
    }
  }
  int num_global = num_checked;
  Mpi::GlobalSum(1, &num_global, comm);
  CHECK((Mpi::Size(comm) == 1 || num_global > 0));
}

TEST_CASE("FaceNbrFieldExchange", "[surfacefunctional][Serial][Parallel][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto elem_type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  auto order = GENERATE(1, 2);
  auto nonconformal = GENERATE(false, true);
  CAPTURE(elem_type, order, nonconformal);

  auto mesh = nonconformal ? MakeNCInterfaceMesh(comm, elem_type)
                           : MakeInterfaceMesh(comm, elem_type);
  auto &pmesh = mesh->Get();

  mfem::ND_FECollection nd_fec(order, pmesh.Dimension());
  mfem::RT_FECollection rt_fec(order - 1, pmesh.Dimension());
  FiniteElementSpace nd_fespace(*mesh, &nd_fec), rt_fespace(*mesh, &rt_fec);

  // Project non-trivial smooth fields. The reference values are computed from the same
  // projected grid functions through the legacy mfem face neighbor dof exchange, so the
  // comparison is exact up to roundoff (not projection error).
  mfem::ParGridFunction E(&nd_fespace.Get()), B(&rt_fespace.Get());
  mfem::VectorFunctionCoefficient fe(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = std::sin(x(1)) + x(2) * x(2);
                                       v(1) = std::cos(x(2)) + x(0);
                                       v(2) = x(0) * x(1) + 1.0;
                                     });
  mfem::VectorFunctionCoefficient fb(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = x(1) * x(2) - 0.5;
                                       v(1) = std::sin(x(0)) - x(2);
                                       v(2) = std::cos(x(1)) + x(0) * x(0);
                                     });
  E.ProjectCoefficient(fe);
  B.ProjectCoefficient(fb);

  // Request E (slot 0) at a few reference points of every ghost element, and both E
  // and B (slot 1) for every other ghost element (exercising the value layouts). The
  // points are valid reference coordinates for both tetrahedra and hexahedra.
  const int num_ghost = pmesh.GetNFaceNeighborElements();
  std::vector<FaceNbrFieldExchange::Request> requests;
  for (int fn = 0; fn < num_ghost; fn++)
  {
    auto &req = requests.emplace_back();
    req.face_nbr_elem = fn;
    req.source_mask = (fn % 2 == 0) ? 0b01u : 0b11u;
    // Deliberately caller-supplied (not a generated trace-map) key. The rebuild below
    // retains it while changing the coordinates, so the exact point signature must
    // distinguish the two process-lifetime fixed rules.
    req.point_key = {17, 29, 4};
    req.pts.resize(4);
    req.pts[0].Set3(0.1, 0.2, 0.3);
    req.pts[1].Set3(0.25, 0.25, 0.25);
    req.pts[2].Set3(0.05, 0.1, 0.7);
    req.pts[3].Set3(0.3, 0.05, 0.05);
  }
  FaceNbrFieldExchange exchange(
      *mesh, {&nd_fespace.Get(), &rt_fespace.Get(), nullptr, nullptr}, requests);
  exchange.Exchange({&E, &B, nullptr, nullptr});

  // Reference: evaluate the ghost elements through the legacy dof exchange.
  E.ExchangeFaceNbrData();
  B.ExchangeFaceNbrData();
  std::vector<double> vals(exchange.Imported().HostRead(),
                           exchange.Imported().HostRead() + exchange.ImportSize());
  mfem::Vector ref(3);
  int num_checked = 0;
  for (std::size_t r = 0; r < requests.size(); r++)
  {
    const auto &req = requests[r];
    for (int s = 0; s < 2; s++)
    {
      const int offset = exchange.ImportOffset(static_cast<int>(r), s);
      if (!(req.source_mask & (1u << s)))
      {
        CHECK(offset < 0);
        continue;
      }
      const auto &U = (s == 0) ? E : B;
      for (std::size_t j = 0; j < req.pts.size(); j++)
      {
        U.GetVectorValue(pmesh.GetNE() + req.face_nbr_elem, req.pts[j], ref);
        for (int c = 0; c < 3; c++)
        {
          CAPTURE(r, s, j, c);
          CHECK(vals[offset + 3 * j + c] == Catch::Approx(ref(c)).margin(1.0e-12));
          num_checked++;
        }
      }
    }
  }
  // With more than one process, the partition must produce at least one ghost element
  // somewhere (the exchange is the point of the test).
  int num_global = num_checked;
  Mpi::GlobalSum(1, &num_global, comm);
  CHECK((Mpi::Size(comm) == 1 || num_global > 0));

  // The field inputs are re-pointed at the sources on each call: scaling the field
  // scales the exchanged values.
  E *= 2.0;
  exchange.Exchange({&E, &B, nullptr, nullptr});
  const double *vals2 = exchange.Imported().HostRead();
  for (std::size_t r = 0; r < requests.size(); r++)
  {
    const int offset = exchange.ImportOffset(static_cast<int>(r), 0);
    for (std::size_t j = 0; j < 3 * requests[r].pts.size(); j++)
    {
      CHECK(vals2[offset + j] == Catch::Approx(2.0 * vals[offset + j]).margin(1.0e-12));
    }
  }

  // Rebuild with the same caller-supplied noncanonical point_key but different points.
  // The exact coordinate-and-weight signature, rather than an exchange ID, must select
  // a distinct process-lifetime IntegrationRule.
  auto requests_2 = requests;
  for (auto &req : requests_2)
  {
    req.pts[0].Set3(0.15, 0.10, 0.20);
    req.pts[1].Set3(0.20, 0.15, 0.30);
    req.pts[2].Set3(0.10, 0.25, 0.20);
    req.pts[3].Set3(0.25, 0.10, 0.15);
  }
  FaceNbrFieldExchange exchange_2(
      *mesh, {&nd_fespace.Get(), &rt_fespace.Get(), nullptr, nullptr}, requests_2);
  exchange_2.Exchange({&E, &B, nullptr, nullptr});
  E.ExchangeFaceNbrData();
  const double *vals_2 = exchange_2.Imported().HostRead();
  for (std::size_t r = 0; r < requests_2.size(); r++)
  {
    const int offset = exchange_2.ImportOffset(static_cast<int>(r), 0);
    for (std::size_t j = 0; j < requests_2[r].pts.size(); j++)
    {
      E.GetVectorValue(pmesh.GetNE() + requests_2[r].face_nbr_elem, requests_2[r].pts[j],
                       ref);
      for (int c = 0; c < 3; c++)
      {
        CHECK(vals_2[offset + 3 * j + c] == Catch::Approx(ref(c)).margin(1.0e-12));
      }
    }
  }
}

}  // namespace palace
