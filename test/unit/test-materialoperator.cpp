// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <vector>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/geodata.hpp"

namespace palace
{
using namespace Catch;

TEST_CASE("MaterialOperator IsIsotropic", "[materialoperator][Serial]")
{
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::TETRAHEDRON));
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh palace_mesh(std::move(par_mesh));

  config::PeriodicBoundaryData periodic;

  SECTION("Trivial isotropic material")
  {
    config::MaterialData material;
    material.attributes = {1};
    // Default values should be isotropic (all eigenvalues = 1).

    MaterialOperator mat_op({material}, periodic, ProblemType::ELECTROSTATIC, palace_mesh);
    REQUIRE(mat_op.IsIsotropic(1) == true);
  }

  SECTION("Non-trivial isotropic material")
  {
    config::MaterialData material;
    material.attributes = {1};
    material.mu_r.s[0] = 2.0;
    material.mu_r.s[1] = 2.0;
    material.mu_r.s[2] = 2.0;

    MaterialOperator mat_op({material}, periodic, ProblemType::ELECTROSTATIC, palace_mesh);
    REQUIRE(mat_op.IsIsotropic(1) == true);
  }

  SECTION("Anisotropic materials")
  {
    config::MaterialData material;
    material.attributes = {1};

    SECTION("Anisotropic permeability")
    {
      material.mu_r.s[0] = 2.0;
      MaterialOperator mat_op({material}, periodic, ProblemType::ELECTROSTATIC,
                              palace_mesh);
      REQUIRE(mat_op.IsIsotropic(1) == false);
    }

    SECTION("Anisotropic permittivity")
    {
      material.epsilon_r.s[1] = 2.0;
      MaterialOperator mat_op({material}, periodic, ProblemType::ELECTROSTATIC,
                              palace_mesh);
      REQUIRE(mat_op.IsIsotropic(1) == false);
    }

    SECTION("Anisotropic loss tangent")
    {
      material.tandelta.s[2] = 0.02;
      MaterialOperator mat_op({material}, periodic, ProblemType::ELECTROSTATIC,
                              palace_mesh);
      REQUIRE(mat_op.IsIsotropic(1) == false);
    }

    SECTION("Anisotropic conductivity")
    {
      material.sigma.s[0] = 1e6;
      material.sigma.s[2] = 2e6;
      MaterialOperator mat_op({material}, periodic, ProblemType::ELECTROSTATIC,
                              palace_mesh);
      REQUIRE(mat_op.IsIsotropic(1) == false);
    }
  }
}

TEST_CASE("MaterialOperator requires materials for retained mesh domains",
          "[materialoperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int size = Mpi::Size(comm);

  // Give every rank two contiguous x-directed coarse cells. Attribute 2 is confined to
  // the first cell on rank 0, with an attribute-1 cell separating it from the partition
  // interface. Thus no other rank sees attribute 2, even through a shared-face neighbor.
  auto serial_mesh = std::make_unique<mfem::Mesh>(mfem::Mesh::MakeCartesian3D(
      2 * size, 1, 1, mfem::Element::TETRAHEDRON, 2.0 * size, 1.0, 1.0));
  std::vector<int> partitioning(serial_mesh->GetNE());
  for (int i = 0; i < serial_mesh->GetNE(); i++)
  {
    const auto *element = serial_mesh->GetElement(i);
    const auto *vertices = element->GetVertices();
    double center_x = 0.0;
    for (int j = 0; j < element->GetNVertices(); j++)
    {
      center_x += serial_mesh->GetVertex(vertices[j])[0];
    }
    center_x /= element->GetNVertices();
    serial_mesh->SetAttribute(i, center_x < 1.0 ? 2 : 1);
    partitioning[i] = static_cast<int>(center_x / 2.0);
  }
  serial_mesh->SetAttributes();
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh, partitioning.data());
  Mesh palace_mesh(std::move(par_mesh));

  const auto &local_attributes = palace_mesh.GetCeedAttributes();
  const bool has_attr2 = local_attributes.find(2) != local_attributes.end();
  CHECK(has_attr2 == (Mpi::Rank(comm) == 0));
  CHECK(palace_mesh.GetNE() > 0);

  config::MaterialData material1;
  material1.attributes = {1};
  config::MaterialData material2;
  material2.attributes = {2};
  config::PeriodicBoundaryData periodic;

  SECTION("All retained attributes have materials")
  {
    CHECK_NOTHROW(MaterialOperator({material1, material2}, periodic,
                                   ProblemType::ELECTROSTATIC, palace_mesh));
  }

  SECTION("Missing material is rejected collectively")
  {
    CHECK_THROWS_WITH(
        MaterialOperator({material1}, periodic, ProblemType::ELECTROSTATIC, palace_mesh),
        Catch::Matchers::ContainsSubstring(
            "Mesh domain attribute 2 has no corresponding entry in "
            "config[\"Domains\"][\"Materials\"]!"));
  }
}

TEST_CASE("MaterialOperator anisotropic wave-number bound", "[materialoperator][Serial]")
{
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::TETRAHEDRON));
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh palace_mesh(std::move(par_mesh));

  config::MaterialData material;
  material.attributes = {1};
  material.mu_r.s = {4.0, 1.0, 1.0};
  material.epsilon_r.s = {1.0, 4.0, 1.0};

  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({material}, periodic, ProblemType::DRIVEN, palace_mesh);

  // Same-axis products max(mu_i * epsilon_i) are only 4, but a wave polarized with E
  // along y and H along x samples mu_x * epsilon_y = 16.
  CHECK(mat_op.GetMaxMuEpsilon() == Approx(16.0));
}

TEST_CASE("MaterialOperator utility functions", "[materialoperator][Serial]")
{
  SECTION("IsOrthonormal")
  {
    config::SymmetricMatrixData<3> data(1.0);

    SECTION("Orthonormal vectors")
    {
      data.v[0] = {1.0, 0.0, 0.0};
      data.v[1] = {0.0, 1.0, 0.0};
      data.v[2] = {0.0, 0.0, 1.0};
      data.s = {1.0, 2.0, 3.0};
      REQUIRE(internal::mat::IsOrthonormal(data));
    }

    SECTION("Non-normalized vectors")
    {
      data.v[0] = {2.0, 0.0, 0.0};
      data.v[1] = {0.0, 1.0, 0.0};
      data.v[2] = {0.0, 0.0, 1.0};
      data.s = {1.0, 1.0, 1.0};
      REQUIRE(!internal::mat::IsOrthonormal(data));
    }

    SECTION("Non-orthogonal vectors")
    {
      data.v[0] = {1.0, 0.0, 0.0};
      data.v[1] = {0.5, 0.866, 0.0};
      data.v[2] = {0.0, 0.0, 1.0};
      data.s = {1.0, 1.0, 1.0};
      REQUIRE(!internal::mat::IsOrthonormal(data));
    }
  }

  SECTION("IsValid")
  {
    config::SymmetricMatrixData<3> data(1.0);

    SECTION("Valid orthonormal matrix with positive eigenvalues")
    {
      data.v[0] = {1.0, 0.0, 0.0};
      data.v[1] = {0.0, 1.0, 0.0};
      data.v[2] = {0.0, 0.0, 1.0};
      data.s = {1.0, 1.0, 1.0};
      REQUIRE(internal::mat::IsValid(data));
    }

    SECTION("Valid orthonormal matrix with different positive eigenvalues")
    {
      data.v[0] = {1.0, 0.0, 0.0};
      data.v[1] = {0.0, 1.0, 0.0};
      data.v[2] = {0.0, 0.0, 1.0};
      data.s = {2.0, 3.0, 4.0};
      REQUIRE(internal::mat::IsValid(data));
    }

    SECTION("Invalid - zero eigenvalue")
    {
      data.v[0] = {1.0, 0.0, 0.0};
      data.v[1] = {0.0, 1.0, 0.0};
      data.v[2] = {0.0, 0.0, 1.0};
      data.s = {1.0, 0.0, 1.0};
      REQUIRE(!internal::mat::IsValid(data));
    }
  }

  SECTION("IsIdentity")
  {
    config::SymmetricMatrixData<3> data(1.0);

    SECTION("Identity matrix")
    {
      data.v[0] = {1.0, 0.0, 0.0};
      data.v[1] = {0.0, 1.0, 0.0};
      data.v[2] = {0.0, 0.0, 1.0};
      data.s = {1.0, 1.0, 1.0};
      REQUIRE(internal::mat::IsIdentity(data));
    }

    SECTION("Non-identity - different eigenvalues")
    {
      data.v[0] = {1.0, 0.0, 0.0};
      data.v[1] = {0.0, 1.0, 0.0};
      data.v[2] = {0.0, 0.0, 1.0};
      data.s = {2.0, 1.0, 1.0};
      REQUIRE(!internal::mat::IsIdentity(data));
    }

    SECTION("Identity but rotated basis")
    {
      data.v[0] = {0.0, 1.0, 0.0};
      data.v[1] = {1.0, 0.0, 0.0};
      data.v[2] = {0.0, 0.0, 1.0};
      data.s = {1.0, 1.0, 1.0};
      REQUIRE(internal::mat::IsIdentity(data));
    }
  }

  SECTION("IsIsotropic")
  {
    config::SymmetricMatrixData<3> data(1.0);

    SECTION("Isotropic material - all eigenvalues equal")
    {
      data.v[0] = {1.0, 0.0, 0.0};
      data.v[1] = {0.0, 1.0, 0.0};
      data.v[2] = {0.0, 0.0, 1.0};
      data.s = {2.5, 2.5, 2.5};
      REQUIRE(internal::mat::IsIsotropic(data));
    }

    SECTION("Anisotropic material - different eigenvalues")
    {
      data.v[0] = {1.0, 0.0, 0.0};
      data.v[1] = {0.0, 1.0, 0.0};
      data.v[2] = {0.0, 0.0, 1.0};
      data.s = {1.0, 2.0, 1.0};
      REQUIRE(!internal::mat::IsIsotropic(data));
    }

    SECTION("Invalid - non-orthonormal but equal eigenvalues")
    {
      data.v[0] = {2.0, 0.0, 0.0};
      data.v[1] = {0.0, 1.0, 0.0};
      data.v[2] = {0.0, 0.0, 1.0};
      data.s = {2.0, 2.0, 2.0};
      REQUIRE(!internal::mat::IsIsotropic(data));
    }
  }
}

TEST_CASE("MaterialOperator ignores ghost-only submesh attributes",
          "[materialoperator][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int rank = Mpi::Rank(comm);
  REQUIRE(Mpi::Size(comm) >= 2);

  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(2, 1, 1, mfem::Element::HEXAHEDRON, 2.0, 1.0, 1.0));
  for (int i = 0; i < serial_mesh->GetNE(); i++)
  {
    serial_mesh->SetAttribute(i, 7);
  }
  serial_mesh->SetAttributes();
  serial_mesh->EnsureNCMesh(true);

  int partitioning[2] = {0, 1};
  auto parent = std::make_unique<mfem::ParMesh>(comm, *serial_mesh, partitioning);
  mfem::Array<mfem::Refinement> refinements;
  if (rank == 0)
  {
    refinements.Append(mfem::Refinement(0));
  }
  parent->GeneralRefinement(refinements, -1);

  mfem::Array<int> surface_attributes;
  surface_attributes.Append(1);
  auto submesh = std::make_unique<mfem::ParSubMesh>(
      mfem::ParSubMesh::CreateFromBoundary(*parent, surface_attributes));
  Mesh palace_mesh(std::move(submesh));
  mesh::RemapSubMeshAttributes(static_cast<mfem::ParSubMesh &>(palace_mesh.Get()));
  palace_mesh.RebuildCeedAttributes();

  const bool owns_elements = rank < 2;
  CHECK((palace_mesh.GetNE() > 0) == owns_elements);
  for (int i = 0; i < palace_mesh.GetNE(); i++)
  {
    CHECK(palace_mesh.Get().GetAttribute(i) == 7);
  }
  const auto &local_attributes = palace_mesh.GetCeedAttributes();
  if (owns_elements)
  {
    CHECK(local_attributes.find(1) != local_attributes.end());
    CHECK(local_attributes.find(7) != local_attributes.end());
    // The remap must refresh the cached distinct-attribute array.
    CHECK(palace_mesh.Get().attributes.Max() == 7);
  }
  else
  {
    CHECK(local_attributes.empty());
  }

  config::MaterialData material;
  material.attributes = {7};
  config::PeriodicBoundaryData periodic;
  CHECK_NOTHROW(MaterialOperator({material}, periodic, ProblemType::DRIVEN, palace_mesh));
}

}  // namespace palace
