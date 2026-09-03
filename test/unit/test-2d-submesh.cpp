// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "fem/interpolator.hpp"
#include "fem/mesh.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

namespace palace
{
using namespace Catch::Matchers;

TEST_CASE("RotateMaterialDefinitions", "[materialoperator][Serial]")
{
  // Create a 2D mesh and anisotropic sapphire-like material, rotate the config-level
  // material definitions into a known local frame, and verify the resulting 2D-native
  // MaterialOperator tensors.
  MPI_Comm comm = Mpi::World();
  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.model.Lc = 1.0;

  auto &material = iodata.domains.materials.emplace_back();
  material.attributes = {1};
  // Anisotropic permittivity: eps = diag(2, 3, 5) in global frame.
  material.epsilon_r.s = {2.0, 3.0, 5.0};

  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(2, 2, mfem::Element::TRIANGLE, false, 1.0, 1.0));
  iodata.NondimensionalizeInputs(serial_mesh);
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  iodata.CheckConfiguration();
  Mesh palace_mesh(std::move(par_mesh));

  {
    // Unrotated: 2D MaterialOperator truncates to leading 2x2 = diag(2, 3); scalar
    // out-of-plane = eps_zz = 5.
    MaterialOperator mat_op(iodata, palace_mesh);
    const auto &eps = mat_op.GetPermittivityReal();
    REQUIRE(eps.SizeI() == 2);
    CHECK_THAT(eps(0, 0, 0), WithinAbs(2.0, 1e-12));
    CHECK_THAT(eps(1, 1, 0), WithinAbs(3.0, 1e-12));
    CHECK_THAT(mat_op.GetPermittivityScalar()(0, 0, 0), WithinAbs(5.0, 1e-12));
  }

  SECTION("Rotation with normal along x-axis")
  {
    // Surface normal = x, tangent frame: e1 = z, e2 = y.
    // Rotated eps_2x2 = R^T diag(2,3,5) R where R = [e1|e2] = [[0,0],[0,1],[1,0]]
    // Result: eps_2x2 = diag(5, 3) (zz, yy components); eps_nn = eps_xx = 2.
    mfem::Vector e1(3), e2(3), normal(3);
    e1 = 0.0;
    e1(2) = 1.0;  // z
    e2 = 0.0;
    e2(1) = 1.0;  // y
    normal = 0.0;
    normal(0) = 1.0;  // x

    RotateMaterialDefinitions(iodata.domains.materials, e1, e2, normal);
    MaterialOperator mat_op(iodata, palace_mesh);

    const auto &eps = mat_op.GetPermittivityReal();
    CHECK_THAT(eps(0, 0, 0), WithinAbs(5.0, 1e-12));  // e1=z component
    CHECK_THAT(eps(1, 1, 0), WithinAbs(3.0, 1e-12));  // e2=y component
    CHECK_THAT(eps(0, 1, 0), WithinAbs(0.0, 1e-12));  // off-diagonal

    // Scalar (out-of-plane) = n^T eps n = eps_xx = 2.0.
    CHECK_THAT(mat_op.GetPermittivityScalar()(0, 0, 0), WithinAbs(2.0, 1e-12));

    // Curl-curl inv permeability (scalar, out-of-plane mu^{-1}).
    // mu = diag(1,1,1) by default, so mu^{-1}_nn = 1.0.
    CHECK_THAT(mat_op.GetCurlCurlInvPermeability()(0, 0, 0), WithinAbs(1.0, 1e-12));
  }

  SECTION("Rotation with normal along z-axis (identity)")
  {
    // Surface normal = z, tangent frame: e1 = x, e2 = y.
    // Rotated eps_2x2 = diag(2, 3) — same as the truncated result.
    mfem::Vector e1(3), e2(3), normal(3);
    e1 = 0.0;
    e1(0) = 1.0;  // x
    e2 = 0.0;
    e2(1) = 1.0;  // y
    normal = 0.0;
    normal(2) = 1.0;  // z

    RotateMaterialDefinitions(iodata.domains.materials, e1, e2, normal);
    MaterialOperator mat_op(iodata, palace_mesh);

    const auto &eps = mat_op.GetPermittivityReal();
    CHECK_THAT(eps(0, 0, 0), WithinAbs(2.0, 1e-12));
    CHECK_THAT(eps(1, 1, 0), WithinAbs(3.0, 1e-12));

    // Scalar (out-of-plane) = eps_zz = 5.0.
    CHECK_THAT(mat_op.GetPermittivityScalar()(0, 0, 0), WithinAbs(5.0, 1e-12));
  }
}

TEST_CASE("Project3Dto2D", "[geodata][Serial]")
{
  // Project known 3D points to 2D using a known tangent frame.
  mfem::Vector centroid(3), e1(3), e2(3);
  centroid = 0.0;
  centroid(0) = 1.0;
  centroid(1) = 2.0;
  centroid(2) = 3.0;

  // Tangent frame: e1 = x, e2 = y (surface normal = z).
  e1 = 0.0;
  e1(0) = 1.0;
  e2 = 0.0;
  e2(1) = 1.0;

  SECTION("Point at centroid projects to origin")
  {
    mfem::Vector p3d(3);
    p3d(0) = 1.0;
    p3d(1) = 2.0;
    p3d(2) = 3.0;
    auto p2d = mesh::Project3Dto2D(p3d, centroid, e1, e2);
    CHECK_THAT(p2d(0), WithinAbs(0.0, 1e-14));
    CHECK_THAT(p2d(1), WithinAbs(0.0, 1e-14));
  }

  SECTION("Offset along e1 gives positive u")
  {
    mfem::Vector p3d(3);
    p3d(0) = 2.5;
    p3d(1) = 2.0;
    p3d(2) = 3.0;  // 1.5 along x from centroid
    auto p2d = mesh::Project3Dto2D(p3d, centroid, e1, e2);
    CHECK_THAT(p2d(0), WithinAbs(1.5, 1e-14));
    CHECK_THAT(p2d(1), WithinAbs(0.0, 1e-14));
  }

  SECTION("Rotated frame: e1 = z, e2 = y (x-normal surface)")
  {
    mfem::Vector re1(3), re2(3);
    re1 = 0.0;
    re1(2) = 1.0;  // z
    re2 = 0.0;
    re2(1) = 1.0;  // y

    mfem::Vector p3d(3);
    p3d(0) = 1.0;
    p3d(1) = 4.0;
    p3d(2) = 5.0;
    auto p2d = mesh::Project3Dto2D(p3d, centroid, re1, re2);
    // u = (p-c)·e1 = (0,2,2)·(0,0,1) = 2
    // v = (p-c)·e2 = (0,2,2)·(0,1,0) = 2
    CHECK_THAT(p2d(0), WithinAbs(2.0, 1e-14));
    CHECK_THAT(p2d(1), WithinAbs(2.0, 1e-14));
  }
}

TEST_CASE("RemapSubMeshBdrAttributes", "[geodata][Serial]")
{
  // Create a 3D mesh, extract a boundary ParSubMesh, remap boundary attributes,
  // and verify the edges get the correct parent boundary face attributes.
  MPI_Comm comm = Mpi::World();

  // MakeCartesian3D: boundary face attributes are 1-6 for the 6 faces of the cube.
  // Face attr assignment: 1=bottom(z=0), 2=front(y=0), 3=right(x=1),
  //                       4=back(y=1), 5=left(x=0), 6=top(z=1)
  // (The exact assignment depends on MFEM's convention.)
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::TETRAHEDRON, 1.0, 1.0, 1.0));
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);

  // Extract submesh from face attr 5 (left face at x=0).
  mfem::Array<int> surface_attrs;
  surface_attrs.Append(5);
  auto par_submesh = mfem::ParSubMesh::CreateFromBoundary(*par_mesh, surface_attrs);

  REQUIRE(par_submesh.GetNE() > 0);

  // Remap boundary attributes: edges shared with other boundary faces should get
  // those faces' attributes, not the surface attribute (5).
  mesh::RemapSubMeshBdrAttributes(par_submesh, surface_attrs);

  // After remapping, NO submesh boundary edge should have the surface attribute (5),
  // because we prefer non-surface attributes.
  bool has_surface_attr = false;
  for (int i = 0; i < par_submesh.GetNBE(); i++)
  {
    if (par_submesh.GetBdrAttribute(i) == 5)
    {
      has_surface_attr = true;
      break;
    }
  }
  CHECK_FALSE(has_surface_attr);

  // All boundary attributes should be from adjacent faces (1, 2, 3, 4, or 6).
  for (int i = 0; i < par_submesh.GetNBE(); i++)
  {
    int attr = par_submesh.GetBdrAttribute(i);
    CHECK(attr != 5);
    CHECK(attr >= 1);
    CHECK(attr <= 6);
  }
}

TEST_CASE("NC parent marker transfer over-constrains port H1 DoFs", "[geodata][Serial]")
{
  auto mesh_path = std::string(PALACE_TEST_DATA_DIR) + "/mesh/nc-boundary-submesh.mesh";
  mfem::Mesh serial_mesh(mesh_path.c_str(), 1, 1, true);
  REQUIRE(serial_mesh.Nonconforming());
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), serial_mesh);

  mfem::Array<int> surface_attrs;
  surface_attrs.Append(14);
  auto submesh = mfem::ParSubMesh::CreateFromBoundary(*par_mesh, surface_attrs);
  REQUIRE(submesh.Nonconforming());
  REQUIRE(submesh.GetNE() > 0);
  mesh::RemapSubMeshBdrAttributes(submesh, surface_attrs);
  mfem::H1_FECollection parent_fec(1, 3), submesh_fec(1, 2);
  mfem::ParFiniteElementSpace parent_fes(par_mesh.get(), &parent_fec);
  mfem::ParFiniteElementSpace submesh_fes(&submesh, &submesh_fec);
  mfem::ParGridFunction parent_marker_field(&parent_fes),
      submesh_marker_field(&submesh_fes);
  auto transfer =
      mfem::ParSubMesh::CreateTransferMap(parent_marker_field, submesh_marker_field);
  mfem::Array<int> attr4, parent_essential, direct_essential;
  attr4.Append(4);
  mfem::Array<int> parent_marker, submesh_marker;
  mesh::AttrToMarker(par_mesh->bdr_attributes.Max(), attr4, parent_marker);
  mesh::AttrToMarker(submesh.bdr_attributes.Max(), attr4, submesh_marker);
  parent_fes.GetEssentialTrueDofs(parent_marker, parent_essential);
  submesh_fes.GetEssentialTrueDofs(submesh_marker, direct_essential);
  mfem::Vector parent_true(parent_fes.GetTrueVSize());
  parent_true = 0.0;
  for (int tdof : parent_essential)
  {
    parent_true[tdof] = 1.0;
  }
  parent_marker_field.SetFromTrueDofs(parent_true);
  transfer.Transfer(parent_marker_field, submesh_marker_field);
  mfem::Vector transferred_true(submesh_fes.GetTrueVSize());
  submesh_marker_field.ParallelProject(transferred_true);
  int transferred_essential = 0;
  for (int i = 0; i < transferred_true.Size(); i++)
  {
    transferred_essential += transferred_true[i] != 0.0;
  }
  // Parent-marker transfer over-constrains two interior H1 DoFs on this NC surface.
  // Essential DoFs must instead be extracted directly from the remapped submesh boundary.
  CHECK(direct_essential.Size() == 6);
  CHECK(transferred_essential == 8);
}

TEST_CASE("Flattened NC port mesh preserves ND DoFs", "[geodata][Serial][Parallel]")
{
  auto mesh_path = std::string(PALACE_TEST_DATA_DIR) + "/mesh/nc-boundary-submesh.mesh";
  mfem::Mesh serial_mesh(mesh_path.c_str(), 1, 1, true);
  auto parent = std::make_unique<mfem::ParMesh>(Mpi::World(), serial_mesh);
  mfem::Array<int> surface_attrs;
  surface_attrs.Append(14);
  auto embedded = mfem::ParSubMesh::CreateFromBoundary(*parent, surface_attrs);
  REQUIRE(embedded.Dimension() == 2);
  REQUIRE(embedded.SpaceDimension() == 3);
  embedded.SetCurvature(3, false, 3, mfem::Ordering::byNODES);
  REQUIRE(embedded.GetNodes()->FESpace()->GetMaxElementOrder() == 3);

  const mfem::Vector normal = mesh::GetSurfaceNormal(embedded);
  mfem::Vector origin(3), e1(3), e2(3);
  origin = 0.0;
  e1 = 0.0;
  e1[1] = 1.0;
  e1.Add(-(e1 * normal), normal);
  e1 /= e1.Norml2();
  e2[0] = normal[1] * e1[2] - normal[2] * e1[1];
  e2[1] = normal[2] * e1[0] - normal[0] * e1[2];
  e2[2] = normal[0] * e1[1] - normal[1] * e1[0];
  e2 /= e2.Norml2();

  mfem::ParMesh flattened(embedded, true);
  mesh::ProjectMeshTo2D(flattened, origin, e1, e2);
  REQUIRE(flattened.SpaceDimension() == 2);

  constexpr int order = 3;
  mfem::ND_FECollection embedded_fec(order, 2), flattened_fec(order, 2);
  mfem::ParFiniteElementSpace embedded_fes(&embedded, &embedded_fec);
  mfem::ParFiniteElementSpace flattened_fes(&flattened, &flattened_fec);
  REQUIRE(embedded_fes.GetVSize() == flattened_fes.GetVSize());
  REQUIRE(embedded_fes.GetTrueVSize() == flattened_fes.GetTrueVSize());

  mfem::VectorFunctionCoefficient embedded_coefficient(
      3,
      [&](const mfem::Vector &, mfem::Vector &value)
      {
        value.SetSize(3);
        for (int d = 0; d < 3; d++)
        {
          value[d] = 1.25 * e1[d] - 0.75 * e2[d];
        }
      });
  mfem::VectorFunctionCoefficient flattened_coefficient(
      2,
      [](const mfem::Vector &, mfem::Vector &value)
      {
        value.SetSize(2);
        value[0] = 1.25;
        value[1] = -0.75;
      });
  mfem::ParGridFunction embedded_field(&embedded_fes), flattened_field(&flattened_fes);
  embedded_field.ProjectCoefficient(embedded_coefficient);
  flattened_field.ProjectCoefficient(flattened_coefficient);
  mfem::ParGridFunction difference(embedded_field);
  difference -= flattened_field;
  double difference_norm = difference.Norml2();
  Mpi::GlobalMax(1, &difference_norm, flattened.GetComm());
  CHECK(difference_norm < 1.0e-12);

#if defined(MFEM_USE_GSLIB)
  mfem::IntegrationPoint ip0, ip1;
  ip0.Set2(0.2, 0.2);
  ip1.Set2(0.6, 0.2);
  mfem::Vector p0(2), p1(2);
  int point_owner = flattened.GetNE() > 0 ? Mpi::Rank(flattened.GetComm())
                                          : Mpi::Size(flattened.GetComm());
  Mpi::GlobalMin(1, &point_owner, flattened.GetComm());
  REQUIRE(point_owner < Mpi::Size(flattened.GetComm()));
  if (Mpi::Rank(flattened.GetComm()) == point_owner)
  {
    auto *transformation = flattened.GetElementTransformation(0);
    transformation->Transform(ip0, p0);
    transformation->Transform(ip1, p1);
  }
  Mpi::Broadcast(2, p0.HostReadWrite(), point_owner, flattened.GetComm());
  Mpi::Broadcast(2, p1.HostReadWrite(), point_owner, flattened.GetComm());
  mfem::FindPointsGSLIB point_locator(flattened.GetComm());
  fem::SetupInterpolator(point_locator, flattened);
  const double integral =
      fem::ComputeLineIntegral(point_locator, p0, p1, flattened_field, 10);
  const double expected = 1.25 * (p1[0] - p0[0]) - 0.75 * (p1[1] - p0[1]);
  CHECK_THAT(integral, WithinAbs(expected, 1.0e-11));
#endif
}

TEST_CASE("Tangent frame from SubMesh extraction", "[geodata][Serial]")
{
  // Build a simple 3D mesh, extract each face with mfem::SubMesh::CreateFromBoundary,
  // project to 2D with mesh::ProjectSubmeshTo2D, and verify the resulting tangent frame
  // (e1, e2, normal) is orthonormal. This mirrors the sequence BoundaryModeSolver::
  // Preprocess runs on the serial mesh before partitioning.

  // Create a unit cube mesh (single hex element).
  mfem::Mesh smesh = mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::HEXAHEDRON);
  smesh.EnsureNodes();

  auto CheckOrthonormal =
      [](const mfem::Vector &e1, const mfem::Vector &e2, const mfem::Vector &n)
  {
    CHECK_THAT(e1.Norml2(), WithinAbs(1.0, 1e-12));
    CHECK_THAT(e2.Norml2(), WithinAbs(1.0, 1e-12));
    CHECK_THAT(e1 * e2, WithinAbs(0.0, 1e-12));
    CHECK_THAT(e1 * n, WithinAbs(0.0, 1e-12));
    CHECK_THAT(e2 * n, WithinAbs(0.0, 1e-12));
  };

  // Unit cube has 6 boundary attributes (one per face). Extract each face as a 2D submesh
  // and verify the tangent frame.
  for (int face_attr = 1; face_attr <= 6; face_attr++)
  {
    SECTION("Face attribute " + std::to_string(face_attr))
    {
      mfem::Array<int> surface_attrs;
      surface_attrs.Append(face_attr);

      auto sub = std::make_unique<mfem::SubMesh>(
          mfem::SubMesh::CreateFromBoundary(smesh, surface_attrs));
      REQUIRE(sub->GetNE() > 0);

      mfem::Vector centroid, e1, e2;
      mfem::Vector normal = mesh::ProjectSubmeshTo2D(*sub, centroid, e1, e2);

      REQUIRE(normal.Size() == 3);
      REQUIRE(e1.Size() == 3);
      REQUIRE(e2.Size() == 3);
      CheckOrthonormal(e1, e2, normal);

      CHECK(sub->Dimension() == 2);
      CHECK(sub->SpaceDimension() == 2);
    }
  }
}

}  // namespace palace
