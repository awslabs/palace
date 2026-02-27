// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"

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

}  // namespace palace
