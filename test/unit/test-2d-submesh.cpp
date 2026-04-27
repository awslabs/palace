// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "fem/mesh.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

namespace palace
{
using namespace Catch::Matchers;

TEST_CASE("RotateMaterialTensors", "[materialoperator][Serial]")
{
  // Create a 2D mesh and MaterialOperator with anisotropic sapphire-like material.
  // Then rotate from global to a known local frame and verify tensor values.
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
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  iodata.NondimensionalizeInputs(*par_mesh);
  iodata.CheckConfiguration();
  Mesh palace_mesh(std::move(par_mesh));
  MaterialOperator mat_op(iodata, palace_mesh);

  // Before rotation: 2D MaterialOperator truncates to leading 2x2 = diag(2, 3).
  const auto &eps_before = mat_op.GetPermittivityReal();
  // Material index 0 (the only material).
  REQUIRE(eps_before.SizeI() == 2);
  CHECK_THAT(eps_before(0, 0, 0), WithinAbs(2.0, 1e-12));
  CHECK_THAT(eps_before(1, 1, 0), WithinAbs(3.0, 1e-12));

  SECTION("Rotation with normal along x-axis")
  {
    // Surface normal = x, tangent frame: e1 = z, e2 = y.
    // Rotated eps_2x2 = R^T diag(2,3,5) R where R = [e1|e2] = [[0,0],[0,1],[1,0]]
    // Result: eps_2x2 = diag(5, 3) (zz, yy components).
    mfem::Vector e1(3), e2(3), normal(3);
    e1 = 0.0;
    e1(2) = 1.0;  // z
    e2 = 0.0;
    e2(1) = 1.0;  // y
    normal = 0.0;
    normal(0) = 1.0;  // x

    mat_op.RotateMaterialTensors(iodata, e1, e2, normal);

    const auto &eps = mat_op.GetPermittivityReal();
    CHECK_THAT(eps(0, 0, 0), WithinAbs(5.0, 1e-12));  // e1=z component
    CHECK_THAT(eps(1, 1, 0), WithinAbs(3.0, 1e-12));  // e2=y component
    CHECK_THAT(eps(0, 1, 0), WithinAbs(0.0, 1e-12));  // off-diagonal

    // Scalar (out-of-plane) = n^T eps n = eps_xx = 2.0
    const auto &eps_s = mat_op.GetPermittivityScalar();
    CHECK_THAT(eps_s(0, 0, 0), WithinAbs(2.0, 1e-12));

    // Curl-curl inv permeability (scalar, out-of-plane mu^{-1}).
    // mu = diag(1,1,1) by default, so mu^{-1}_nn = 1.0.
    const auto &muinv_s = mat_op.GetCurlCurlInvPermeability();
    CHECK_THAT(muinv_s(0, 0, 0), WithinAbs(1.0, 1e-12));
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

    mat_op.RotateMaterialTensors(iodata, e1, e2, normal);

    const auto &eps = mat_op.GetPermittivityReal();
    CHECK_THAT(eps(0, 0, 0), WithinAbs(2.0, 1e-12));
    CHECK_THAT(eps(1, 1, 0), WithinAbs(3.0, 1e-12));

    // Scalar (out-of-plane) = eps_zz = 5.0
    const auto &eps_s = mat_op.GetPermittivityScalar();
    CHECK_THAT(eps_s(0, 0, 0), WithinAbs(5.0, 1e-12));
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

TEST_CASE("Tangent frame from SubMesh extraction", "[geodata][Serial]")
{
  // Build a simple 3D mesh, extract each face with mfem::SubMesh::CreateFromBoundary,
  // project to 2D with mesh::ProjectSubmeshTo2D, and verify the resulting tangent frame
  // (e1, e2, normal) is orthonormal. This mirrors the sequence BoundaryModeSolver::
  // PreprocessMesh runs on the serial mesh before partitioning.

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
      mfem::Vector normal = mesh::ProjectSubmeshTo2D(*sub, &centroid, &e1, &e2);

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
