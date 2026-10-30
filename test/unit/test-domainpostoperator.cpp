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
#include "fem/integrator.hpp"
#include "models/domainpostoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/postoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/constants.hpp"
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

  // Initialize quadrature order for tests using granular constructors.
  constexpr int order = 1;
  fem::DefaultIntegrationOrder::p_trial = order;

  Units units(0.496, 1.453);

  config::MaterialData material;
  material.attributes = {1};
  std::vector<config::MaterialData> materials = {material};

  config::PeriodicBoundaryData periodic;
  config::DomainPostData postpro;

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

  mfem::ND_FECollection nd_fec(1, dim);
  FiniteElementSpace nd_fespace(palace_mesh, &nd_fec);
  mfem::RT_FECollection rt_fec(1, dim);
  FiniteElementSpace rt_fespace(palace_mesh, &rt_fec);
  MaterialOperator mat_op(materials, periodic, ProblemType::ELECTROSTATIC, palace_mesh);

  DomainPostOperator dom_post_op(postpro, mat_op, nd_fespace, rt_fespace);

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

TEST_CASE("DomainPostOperator - PML Regions", "[domainpostoperator][pml][Serial][Parallel]")
{
  // Box [0, 1] x [0, 1] x [0, 2] with a PML region z > 1 (attribute 2, εᵣ = 4) on top of
  // the physical region (attribute 1, εᵣ = 1), and a uniform field. The total energy
  // excludes the PML region, whose energies (of the domain postprocessing and of the PML
  // regions) are computed with the background material of the layer.
  constexpr int order = 1;
  fem::DefaultIntegrationOrder::p_trial = order;
  MPI_Comm comm = Mpi::World();
  mfem::Mesh serial_mesh =
      mfem::Mesh::MakeCartesian3D(2, 2, 4, mfem::Element::HEXAHEDRON, 1.0, 1.0, 2.0);
  for (int i = 0; i < serial_mesh.GetNE(); i++)
  {
    mfem::Vector center;
    serial_mesh.GetElementCenter(i, center);
    serial_mesh.SetAttribute(i, (center(2) > 1.0) ? 2 : 1);
  }
  serial_mesh.SetAttributes();
  Mesh palace_mesh(std::make_unique<mfem::ParMesh>(comm, serial_mesh));

  config::MaterialData vacuum, substrate;
  vacuum.attributes = {1};
  substrate.attributes = {2};
  substrate.epsilon_r = 4.0;
  config::PMLData pml;
  pml.attributes = {2};
  pml.autodetect_geometry = false;
  pml.directions = {false, false, false, false, false, true};
  pml.thickness = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
  pml.reference_frequency = 1.0;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({vacuum, substrate}, periodic, ProblemType::DRIVEN, palace_mesh,
                          {pml}, order);
  REQUIRE(mat_op.HasPML());

  config::DomainPostData postpro;
  config::DomainEnergyData energy_1, energy_2;
  energy_1.attributes = {1};
  energy_2.attributes = {2};
  postpro.energy = {{1, energy_1}, {2, energy_2}};

  mfem::ND_FECollection nd_fec(order, 3);
  FiniteElementSpace nd_fespace(palace_mesh, &nd_fec);
  mfem::RT_FECollection rt_fec(order - 1, 3);
  FiniteElementSpace rt_fespace(palace_mesh, &rt_fec);
  DomainPostOperator dom_post_op(postpro, mat_op, nd_fespace, rt_fespace);
  REQUIRE(dom_post_op.HasPML());

  // Uniform fields E = (1, 2, 0) (complex, with E_imag = E_real / 2) and B = (0, 0, 3).
  auto E = std::make_unique<GridFunction>(nd_fespace, true);
  mfem::Vector e_vec({1.0, 2.0, 0.0});
  mfem::VectorConstantCoefficient e_coeff(e_vec);
  E->Real().ProjectCoefficient(e_coeff);
  E->Imag() = E->Real();
  E->Imag() *= 0.5;
  auto B = std::make_unique<GridFunction>(rt_fespace, false);
  mfem::Vector b_vec({0.0, 0.0, 3.0});
  mfem::VectorConstantCoefficient b_coeff(b_vec);
  B->Real().ProjectCoefficient(b_coeff);

  // E_elec = 1/2 εᵣ |E|² V and E_mag = 1/2 |B|² V (unit volumes, μᵣ = 1).
  const double e2 = 1.25 * 5.0, b2 = 9.0;
  CHECK_THAT(dom_post_op.GetElectricFieldEnergy(*E), WithinRel(0.5 * e2, 1.0e-12));
  CHECK_THAT(dom_post_op.GetDomainElectricFieldEnergy(1, *E), WithinRel(0.5 * e2, 1.0e-12));
  CHECK_THAT(dom_post_op.GetDomainElectricFieldEnergy(2, *E),
             WithinRel(0.5 * 4.0 * e2, 1.0e-12));
  CHECK_THAT(dom_post_op.GetPMLElectricFieldEnergy(*E), WithinRel(0.5 * 4.0 * e2, 1.0e-12));
  CHECK_THAT(dom_post_op.GetMagneticFieldEnergy(*B), WithinRel(0.5 * b2, 1.0e-12));
  CHECK_THAT(dom_post_op.GetDomainMagneticFieldEnergy(2, *B), WithinRel(0.5 * b2, 1.0e-12));
  CHECK_THAT(dom_post_op.GetPMLMagneticFieldEnergy(*B), WithinRel(0.5 * b2, 1.0e-12));
}
}  // namespace palace
