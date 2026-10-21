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
#include "utils/permittivitymodel.hpp"
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
  double energy_nondim = dom_post_op.GetElectricFieldEnergy(*E_field, 0.0);
  double energy_SI = units.Dimensionalize<Units::ValueType::ENERGY>(energy_nondim);

  // Expected: U = (ε₀ * E₀² * sx * sy * sz) / 2.
  double sx_SI = units.Dimensionalize<Units::ValueType::LENGTH>(sx);
  double sy_SI = units.Dimensionalize<Units::ValueType::LENGTH>(sy);
  double sz_SI = units.Dimensionalize<Units::ValueType::LENGTH>(sz);
  double expected_energy_SI =
      0.5 * electromagnetics::epsilon0_ * E0_SI * E0_SI * sx_SI * sy_SI * sz_SI;

  CHECK_THAT(energy_SI, WithinRel(expected_energy_SI, 0.01));
}

TEST_CASE("DomainPostOperator - Frequency-dependent electric energy",
          "[domainpostoperator][Serial][Parallel]")
{
  // A frequency-dependent material stores the electric energy of a constant permittivity
  // ε∞ + w(ω), with w the stored-energy susceptibility of its model, over the entire domain
  // and over a postprocessing domain that covers it in part. Its field energy at the
  // material's permittivity is that of the constant ε∞ + Re χ(iω).
  constexpr int order = 1;
  fem::DefaultIntegrationOrder::p_trial = order;
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(4, 2, 2, mfem::Element::HEXAHEDRON, 2.0, 1.0, 1.0));
  for (int i = 0; i < serial_mesh->GetNE(); i++)
  {
    mfem::Vector center;
    serial_mesh->GetElementCenter(i, center);
    serial_mesh->SetAttribute(i, 1 + static_cast<int>(2.0 * center(0)));  // 1 to 4
  }
  serial_mesh->SetAttributes();
  Mesh palace_mesh(std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh));
  mfem::ND_FECollection nd_fec(order, 3);
  FiniteElementSpace nd_fespace(palace_mesh, &nd_fec);
  mfem::RT_FECollection rt_fec(order - 1, 3);
  FiniteElementSpace rt_fespace(palace_mesh, &rt_fec);

  // Attributes 1 and 2 are a Drude and Lorentz material, with ω between their rates, and 3
  // and 4 a constant one. The postprocessing domain is attributes 2 and 3.
  const double omega = 1.7;
  config::MaterialData dispersive, constant;
  dispersive.attributes = {1, 2};
  dispersive.epsilon_r = config::SymmetricMatrixData<3>(2.0);
  dispersive.permittivity_terms = {config::DrudeData{3.0, 0.2},
                                   config::LorentzData{0.8, 2.5, 0.3}};
  constant.attributes = {3, 4};
  constant.epsilon_r = config::SymmetricMatrixData<3>(3.0);
  auto reference = dispersive;
  reference.permittivity_terms.clear();
  using ElectricEnergy = DomainPostOperator::ElectricEnergy;
  const auto type = GENERATE(ElectricEnergy::STORED, ElectricEnergy::MATERIAL);
  CAPTURE(static_cast<int>(type));
  const double w =
      (type == ElectricEnergy::STORED)
          ? permittivity::StoredEnergySusceptibility(dispersive.permittivity_terms, omega)
          : permittivity::Susceptibility(dispersive.permittivity_terms, {0.0, omega})
                .real();
  REQUIRE(w != 0.0);
  reference.epsilon_r = config::SymmetricMatrixData<3>(2.0 + w);
  config::DomainPostData postpro;
  postpro.energy[1].attributes = {2, 3};
  config::PeriodicBoundaryData periodic;
  MaterialOperator dispersive_op({dispersive, constant}, periodic, ProblemType::DRIVEN,
                                 palace_mesh);
  MaterialOperator reference_op({reference, constant}, periodic, ProblemType::DRIVEN,
                                palace_mesh);
  DomainPostOperator dispersive_post(postpro, dispersive_op, nd_fespace, rt_fespace);
  DomainPostOperator reference_post(postpro, reference_op, nd_fespace, rt_fespace);

  GridFunction E(nd_fespace, true);
  mfem::VectorFunctionCoefficient Er(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = 1.0 + x(1);
                                       v(1) = x(0) * x(2);
                                       v(2) = std::sin(x(0));
                                     }),
      Ei(3,
         [](const mfem::Vector &x, mfem::Vector &v)
         {
           v(0) = x(2);
           v(1) = 0.5;
           v(2) = x(0) - x(1);
         });
  E.Real().ProjectCoefficient(Er);
  E.Imag().ProjectCoefficient(Ei);
  const double total = reference_post.GetElectricFieldEnergy(E, omega);
  const double domain = reference_post.GetDomainElectricFieldEnergy(1, E, omega);
  CHECK_THAT(dispersive_post.GetElectricFieldEnergy(E, omega, type),
             WithinRel(total, 1.0e-12));
  CHECK_THAT(dispersive_post.GetDomainElectricFieldEnergy(1, E, omega, type),
             WithinRel(domain, 1.0e-12));
  // The ε∞ energy alone differs.
  auto bare = dispersive;
  bare.permittivity_terms.clear();
  MaterialOperator bare_op({bare, constant}, periodic, ProblemType::DRIVEN, palace_mesh);
  DomainPostOperator bare_post(postpro, bare_op, nd_fespace, rt_fespace);
  CHECK(std::abs(bare_post.GetElectricFieldEnergy(E, omega, type) - total) > 0.1 * total);
}
}  // namespace palace
