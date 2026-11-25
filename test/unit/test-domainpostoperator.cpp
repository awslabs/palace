// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <numbers>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <mfem/fem/coefficient.hpp>
#include <mfem/linalg/vector.hpp>
#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "linalg/vector.hpp"
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

  // The time-averaged E_elec = 1/4 εᵣ |E|² V of the complex phasor E, and the energy
  // E_mag = 1/2 |B|² V of the real field B (unit volumes, μᵣ = 1). The materials do not
  // depend on frequency.
  const double e2 = 1.25 * 5.0, b2 = 9.0, omega = 1.0;
  CHECK_THAT(dom_post_op.GetElectricFieldEnergy(*E, omega), WithinRel(0.25 * e2, 1.0e-12));
  CHECK_THAT(dom_post_op.GetDomainElectricFieldEnergy(1, *E, omega),
             WithinRel(0.25 * e2, 1.0e-12));
  CHECK_THAT(dom_post_op.GetDomainElectricFieldEnergy(2, *E, omega),
             WithinRel(0.25 * 4.0 * e2, 1.0e-12));
  CHECK_THAT(dom_post_op.GetPMLElectricFieldEnergy(*E),
             WithinRel(0.25 * 4.0 * e2, 1.0e-12));
  CHECK_THAT(dom_post_op.GetMagneticFieldEnergy(*B), WithinRel(0.5 * b2, 1.0e-12));
  CHECK_THAT(dom_post_op.GetDomainMagneticFieldEnergy(2, *B), WithinRel(0.5 * b2, 1.0e-12));
  CHECK_THAT(dom_post_op.GetPMLMagneticFieldEnergy(*B), WithinRel(0.5 * b2, 1.0e-12));
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

TEST_CASE("DomainPostOperator - 2D Magnetic Energy",
          "[domainpostoperator][Serial][Parallel]")
{
  // In 2D, B = curl E is a scalar in an L2 space with INTEGRAL map type, B = B̂ / |J|. The
  // time-averaged magnetic energy 1/4 ∫ μ⁻¹ |B|² dA of the complex phasor B must not depend
  // on the element Jacobians. The rectangle [0, 1] x [0, 0.5] has μ_r = 1 for x < 1/2
  // (attribute 1, also the energy region) and μ_r = 2 for x > 1/2 (attribute 2).
  const int n = GENERATE(2, 4);
  const bool distort = GENERATE(false, true);
  CAPTURE(n, distort);
  fem::DefaultIntegrationOrder::p_trial = 2;
  auto serial = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(n, n, mfem::Element::QUADRILATERAL, false, 1.0, 0.5));
  for (int i = 0; i < serial->GetNE(); i++)
  {
    mfem::Vector center;
    serial->GetElementCenter(i, center);
    serial->SetAttribute(i, (center(0) < 0.5) ? 1 : 2);
  }
  serial->SetAttributes();
  if (distort)
  {
    // Move the interior vertices, so the quadrilaterals are not parallelograms. The
    // boundary vertices do not move, so the area is unchanged.
    serial->Transform(
        [](const mfem::Vector &x, mfem::Vector &y)
        {
          y = x;
          y(0) += 0.1 * std::sin(std::numbers::pi * x(0)) *
                  std::sin(2.0 * std::numbers::pi * x(1));
          y(1) += 0.05 * std::sin(2.0 * std::numbers::pi * x(0)) *
                  std::sin(2.0 * std::numbers::pi * x(1));
        });
  }
  if (Mpi::Size(Mpi::World()) > serial->GetNE())
  {
    SKIP("Mesh has fewer elements than processes");
  }
  Mesh mesh(std::make_unique<mfem::ParMesh>(Mpi::World(), *serial));
  mfem::ND_FECollection nd_fec(2, 2);
  mfem::L2_FECollection l2_fec(1, 2, mfem::BasisType::GaussLegendre,
                               mfem::FiniteElement::INTEGRAL);
  FiniteElementSpace nd_fespace(mesh, &nd_fec), l2_fespace(mesh, &l2_fec);
  std::vector<config::MaterialData> materials(2);
  materials[0].attributes = {1};
  materials[1].attributes = {2};
  materials[1].mu_r = 2.0;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op(materials, periodic, ProblemType::EIGENMODE, mesh);
  config::DomainPostData postpro;
  postpro.energy[1].attributes = {1};
  DomainPostOperator dom_post_op(postpro, mat_op, nd_fespace, l2_fespace);

  if (!distort)
  {
    // Project the physical field B = 1 + i x, which the space represents exactly on
    // parallelograms. Region 1: 1/4 ∫ (1 + x²) dA = 1/16 + 1/192. Region 2: 1/8 ∫ (1 + x²)
    // dA = 1/32 + 7/384.
    GridFunction B(l2_fespace, true);
    mfem::ConstantCoefficient one(1.0);
    mfem::FunctionCoefficient x([](const mfem::Vector &x) { return x(0); });
    B.Real().ProjectCoefficient(one);
    B.Imag().ProjectCoefficient(x);
    const double region = 1.0 / 16.0 + 1.0 / 192.0;
    const double total = region + 1.0 / 32.0 + 7.0 / 384.0;
    CHECK_THAT(dom_post_op.GetMagneticFieldEnergy(B), WithinRel(total, 1.0e-12));
    CHECK_THAT(dom_post_op.GetDomainMagneticFieldEnergy(1, B), WithinRel(region, 1.0e-12));
  }

  // The energy of B = curl E must match the time-averaged curl-curl quadratic 1/4 Eᴴ K E of
  // the system matrix on any mesh.
  GridFunction E(nd_fespace, true), B(l2_fespace, true);
  mfem::VectorFunctionCoefficient Er(2,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = -x(1) + x(1) * x(1);
                                       v(1) = x(0) * x(0);
                                     });
  mfem::VectorFunctionCoefficient Ei(2,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = std::sin(std::numbers::pi * x(1));
                                       v(1) =
                                           std::cos(std::numbers::pi * x(0)) + x(0) * x(1);
                                     });
  E.Real().ProjectCoefficient(Er);
  E.Imag().ProjectCoefficient(Ei);
  mfem::DiscreteLinearOperator curl(&nd_fespace.Get(), &l2_fespace.Get());
  curl.AddDomainInterpolator(new mfem::CurlInterpolator);
  curl.Assemble();
  curl.Finalize();
  curl.Mult(E.Real(), B.Real());
  curl.Mult(E.Imag(), B.Imag());
  auto CurlCurlEnergy = [&](bool restrict)
  {
    MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                           mat_op.GetCurlCurlInvPermeability());
    if (restrict)
    {
      muinv_func.RestrictCoefficient(
          mat_op.GetCeedAttributes(postpro.energy[1].attributes));
    }
    BilinearForm k(nd_fespace);
    k.AddDomainIntegrator<CurlCurlIntegrator>(muinv_func);
    auto K = k.PartialAssemble();
    Vector KE(K->Height());
    KE.UseDevice(true);
    K->Mult(E.Real(), KE);
    double dot = linalg::LocalDot(E.Real(), KE);
    K->Mult(E.Imag(), KE);
    dot += linalg::LocalDot(E.Imag(), KE);
    Mpi::GlobalSum(1, &dot, Mpi::World());
    return 0.25 * dot;
  };
  const double total = CurlCurlEnergy(false), region = CurlCurlEnergy(true);
  REQUIRE(region > 0.0);
  REQUIRE(total > region);
  CHECK_THAT(dom_post_op.GetMagneticFieldEnergy(B), WithinRel(total, 1.0e-12));
  CHECK_THAT(dom_post_op.GetDomainMagneticFieldEnergy(1, B), WithinRel(region, 1.0e-12));
}

TEST_CASE("DomainPostOperator - Mode analysis energies",
          "[domainpostoperator][Serial][Parallel]")
{
  // The energies of a boundary mode include the out-of-plane E and the in-plane B. Uniform
  // fields on the rectangle [0, 1] x [0, 0.5]: attribute 1 (x < 1/2, also the energy
  // region) is a frequency-dependent material with ε∞ = 2, and attribute 2 a constant one
  // with ε_r = diag(3, 3, 5) and μ_r = diag(2, 2, 4), whose out-of-plane components weight
  // En and Bz. E = (1, 2, 3) + i (0.5, 0, -1), |Et|² = 5.25 and |En|² = 10, and
  // z_hat x Bt = (0.3, 0.4) and Bz = 0.7, |Bt|² = 0.25 and |Bz|² = 0.49; each region has
  // area 1/4.
  constexpr int order = 1;
  fem::DefaultIntegrationOrder::p_trial = order;
  auto serial = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(4, 4, mfem::Element::QUADRILATERAL, false, 1.0, 0.5));
  for (int i = 0; i < serial->GetNE(); i++)
  {
    mfem::Vector center;
    serial->GetElementCenter(i, center);
    serial->SetAttribute(i, (center(0) < 0.5) ? 1 : 2);
  }
  serial->SetAttributes();
  Mesh mesh(std::make_unique<mfem::ParMesh>(Mpi::World(), *serial));
  mfem::ND_FECollection nd_fec(order, 2);
  mfem::H1_FECollection h1_fec(order, 2);
  mfem::L2_FECollection l2_fec(order - 1, 2, mfem::BasisType::GaussLegendre,
                               mfem::FiniteElement::INTEGRAL);
  FiniteElementSpace nd_fespace(mesh, &nd_fec), h1_fespace(mesh, &h1_fec),
      l2_fespace(mesh, &l2_fec);
  std::vector<config::MaterialData> materials(2);
  materials[0].attributes = {1};
  materials[0].epsilon_r = config::SymmetricMatrixData<3>(2.0);
  materials[0].permittivity_terms = {config::DrudeData{3.0, 0.2},
                                     config::LorentzData{0.8, 2.5, 0.3}};
  materials[1].attributes = {2};
  materials[1].epsilon_r.s = {3.0, 3.0, 5.0};
  materials[1].mu_r.s = {2.0, 2.0, 4.0};
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op(materials, periodic, ProblemType::BOUNDARYMODE, mesh);
  config::DomainPostData postpro;
  postpro.energy[1].attributes = {1};
  DomainPostOperator dom_post_op(postpro, mat_op, nd_fespace, h1_fespace, l2_fespace);

  GridFunction Et(nd_fespace, true), En(h1_fespace, true), zxBt(nd_fespace, true),
      Bz(l2_fespace, true);
  mfem::VectorConstantCoefficient Et_r(mfem::Vector({1.0, 2.0})),
      Et_i(mfem::Vector({0.5, 0.0})), zxBt_r(mfem::Vector({0.3, 0.4})),
      zero_t(mfem::Vector({0.0, 0.0}));
  mfem::ConstantCoefficient En_r(3.0), En_i(-1.0), Bz_r(0.7), zero(0.0);
  Et.Real().ProjectCoefficient(Et_r);
  Et.Imag().ProjectCoefficient(Et_i);
  En.Real().ProjectCoefficient(En_r);
  En.Imag().ProjectCoefficient(En_i);
  zxBt.Real().ProjectCoefficient(zxBt_r);
  zxBt.Imag().ProjectCoefficient(zero_t);
  Bz.Real().ProjectCoefficient(Bz_r);
  Bz.Imag().ProjectCoefficient(zero);

  const double omega = 1.7;
  const double w_stored =
      permittivity::StoredEnergySusceptibility(materials[0].permittivity_terms, omega);
  const double w_material =
      permittivity::Susceptibility(materials[0].permittivity_terms, {0.0, omega}).real();
  const double E2 = 0.25 * 0.25 * 15.25, B2 = 0.25 * 0.25 * 0.74;
  const double E2_aniso = 0.25 * 0.25 * (5.25 * 3.0 + 10.0 * 5.0),
               B2_aniso = 0.25 * 0.25 * (0.25 / 2.0 + 0.49 / 4.0);
  const auto [stored, material] = dom_post_op.GetElectricFieldEnergies(Et, omega, &En);
  REQUIRE(stored);
  CHECK_THAT(*stored, WithinRel(E2 * (2.0 + w_stored) + E2_aniso, 1.0e-12));
  CHECK_THAT(material, WithinRel(E2 * (2.0 + w_material) + E2_aniso, 1.0e-12));
  CHECK_THAT(dom_post_op.GetHighFrequencyElectricFieldEnergy(Et, &En),
             WithinRel(E2 * 2.0 + E2_aniso, 1.0e-12));
  using ElectricEnergy = DomainPostOperator::ElectricEnergy;
  CHECK_THAT(
      dom_post_op.GetDomainElectricFieldEnergy(1, Et, omega, ElectricEnergy::STORED, &En),
      WithinRel(E2 * (2.0 + w_stored), 1.0e-12));
  CHECK_THAT(
      dom_post_op.GetDomainElectricFieldEnergy(1, Et, omega, ElectricEnergy::MATERIAL, &En),
      WithinRel(E2 * (2.0 + w_material), 1.0e-12));
  CHECK_THAT(dom_post_op.GetMagneticFieldEnergy(Bz, &zxBt),
             WithinRel(B2 + B2_aniso, 1.0e-12));
  CHECK_THAT(dom_post_op.GetDomainMagneticFieldEnergy(1, Bz, &zxBt),
             WithinRel(B2, 1.0e-12));
}

}  // namespace palace
