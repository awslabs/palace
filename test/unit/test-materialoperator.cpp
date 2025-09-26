// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

namespace palace
{
using namespace Catch;

TEST_CASE("MaterialOperator IsIsotropic", "[materialoperator][Serial]")
{
  Units units(1.0, 1.0);

  auto serial_mesh = std::make_unique<mfem::Mesh>(
      std::string(PALACE_TEST_MESH_DIR "/gmsh/sphere.msh"), 1, 1);
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh palace_mesh(std::move(par_mesh));

  SECTION("Trivial isotropic material")
  {
    IoData iodata(units);
    auto &material = iodata.domains.materials.emplace_back();
    material.attributes = {1};
    // Default values should be isotropic (all eigenvalues = 1).

    MaterialOperator mat_op(iodata, palace_mesh);
    REQUIRE(mat_op.IsIsotropic(1) == true);
  }

  SECTION("Non-trivial isotropic material")
  {
    IoData iodata(units);
    auto &material = iodata.domains.materials.emplace_back();
    material.attributes = {1};
    material.mu_r.s[0] = 2.0;
    material.mu_r.s[1] = 2.0;
    material.mu_r.s[2] = 2.0;

    MaterialOperator mat_op(iodata, palace_mesh);
    REQUIRE(mat_op.IsIsotropic(1) == true);
  }

  SECTION("Anisotropic materials")
  {
    IoData iodata(units);
    auto &material = iodata.domains.materials.emplace_back();
    material.attributes = {1};

    SECTION("Anisotropic permeability")
    {
      material.mu_r.s[0] = 2.0;
      MaterialOperator mat_op(iodata, palace_mesh);
      REQUIRE(mat_op.IsIsotropic(1) == false);
    }

    SECTION("Anisotropic permittivity")
    {
      material.epsilon_r.s[1] = 2.0;
      MaterialOperator mat_op(iodata, palace_mesh);
      REQUIRE(mat_op.IsIsotropic(1) == false);
    }

    SECTION("Anisotropic loss tangent")
    {
      material.tandelta.s[2] = 0.02;
      MaterialOperator mat_op(iodata, palace_mesh);
      REQUIRE(mat_op.IsIsotropic(1) == false);
    }

    SECTION("Anisotropic conductivity")
    {
      material.sigma.s[0] = 1e6;
      material.sigma.s[2] = 2e6;
      MaterialOperator mat_op(iodata, palace_mesh);
      REQUIRE(mat_op.IsIsotropic(1) == false);
    }
  }
}

}  // namespace palace
