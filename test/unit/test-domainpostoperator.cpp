// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <mfem/fem/coefficient.hpp>
#include <mfem/linalg/vector.hpp>
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "models/domainpostoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/postoperator.hpp"
#include "utils/communication.hpp"
#include "utils/constants.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

namespace palace
{
using namespace Catch;
using namespace Catch::Matchers;

TEST_CASE("DomainPostOperator - Electric Energy Units", "[domainpostoperator][Serial]")
{
  // Checking units by comparing the total electric energy in a cube with a
  // uniform field along the z axis.
  //
  // TODO: This test can be expanded/improved to be a more robust test for the
  // actual function, not just the units.

  Units units(0.496, 1.453);
  IoData iodata = IoData(units);
  iodata.model.Lc = 1.0;  // Keep spatial conversions simple.

  auto &material = iodata.domains.materials.emplace_back();
  material.attributes = {1};

  // Make mesh for a cube [0, sz] x [0, sy] x [0, sz].
  constexpr double sx = 1.1, sy = 2.5, sz = 3.8;
  MPI_Comm comm = Mpi::World();
  constexpr int resolution = 20;
  std::unique_ptr<mfem::Mesh> serial_mesh =
      std::make_unique<mfem::Mesh>(mfem::Mesh::MakeCartesian3D(
          resolution, resolution, resolution, mfem::Element::TETRAHEDRON, sx, sy, sz));
  const int dim = serial_mesh->Dimension();
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  Mesh palace_mesh(std::move(par_mesh));

  iodata.NondimensionalizeInputs(palace_mesh);
  iodata.CheckConfiguration();  // initializes quadrature

  mfem::ND_FECollection nd_fec(1, dim);
  FiniteElementSpace nd_fespace(palace_mesh, &nd_fec);
  mfem::RT_FECollection rt_fec(1, dim);
  FiniteElementSpace rt_fespace(palace_mesh, &rt_fec);
  MaterialOperator mat_op(iodata, palace_mesh);

  DomainPostOperator dom_post_op(iodata, mat_op, nd_fespace, rt_fespace);

  // Create uniform electric field E = E0 * z_hat (SI units: V/m).
  constexpr double E0_SI = 1.0e6;  // 1 MV/m

  auto E_uniform = [&](const mfem::Vector &x, mfem::Vector &E)
  {
    E(0) = 0.0;
    E(1) = 0.0;
    E(2) = units.Nondimensionalize<Units::ValueType::FIELD_E>(E0_SI);
  };

  constexpr bool complex = false;
  auto E_field = std::make_unique<GridFunction>(nd_fespace, complex);
  mfem::VectorFunctionCoefficient E_coeff(dim, E_uniform);
  E_field->Real().ProjectCoefficient(E_coeff);

  // Compute electric field energy.
  double energy_nondim = dom_post_op.GetElectricFieldEnergy(*E_field);
  double energy_SI = units.Dimensionalize<Units::ValueType::ENERGY>(energy_nondim);

  // Expected: U = (ε₀ * E₀² * sx * sy * sz) / 2.
  double sx_SI = units.Dimensionalize<Units::ValueType::LENGTH>(sx);
  double sy_SI = units.Dimensionalize<Units::ValueType::LENGTH>(sy);
  double sz_SI = units.Dimensionalize<Units::ValueType::LENGTH>(sz);
  double expected_energy_SI =
      0.5 * electromagnetics::epsilon0_ * E0_SI * E0_SI * sx_SI * sy_SI * sz_SI;

  CHECK_THAT(energy_SI, WithinRel(expected_energy_SI, 0.01));
}
}  // namespace palace
