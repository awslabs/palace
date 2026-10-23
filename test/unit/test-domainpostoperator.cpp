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

TEST_CASE("DomainPostOperator - 2D Magnetic Energy",
          "[domainpostoperator][Serial][Parallel]")
{
  // In 2D, B = curl E is a scalar in an L2 space with INTEGRAL map type, B = B̂ / |J|. The
  // magnetic energy 1/2 ∫ μ⁻¹ |B|² dA must not depend on the element Jacobians. The
  // rectangle [0, 1] x [0, 0.5] has μ_r = 1 for x < 1/2 (attribute 1, also the energy
  // region) and μ_r = 2 for x > 1/2 (attribute 2).
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
    // parallelograms. Region 1: 1/2 ∫ (1 + x²) dA = 1/8 + 1/96. Region 2: 1/4 ∫ (1 + x²) dA
    // = 1/16 + 7/192.
    GridFunction B(l2_fespace, true);
    mfem::ConstantCoefficient one(1.0);
    mfem::FunctionCoefficient x([](const mfem::Vector &x) { return x(0); });
    B.Real().ProjectCoefficient(one);
    B.Imag().ProjectCoefficient(x);
    const double region = 1.0 / 8.0 + 1.0 / 96.0;
    const double total = region + 1.0 / 16.0 + 7.0 / 192.0;
    CHECK_THAT(dom_post_op.GetMagneticFieldEnergy(B), WithinRel(total, 1.0e-12));
    CHECK_THAT(dom_post_op.GetDomainMagneticFieldEnergy(1, B), WithinRel(region, 1.0e-12));
  }

  // The energy of B = curl E must match the curl-curl quadratic 1/2 Eᴴ K E of the system
  // matrix on any mesh.
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
    return 0.5 * dot;
  };
  const double total = CurlCurlEnergy(false), region = CurlCurlEnergy(true);
  REQUIRE(region > 0.0);
  REQUIRE(total > region);
  CHECK_THAT(dom_post_op.GetMagneticFieldEnergy(B), WithinRel(total, 1.0e-12));
  CHECK_THAT(dom_post_op.GetDomainMagneticFieldEnergy(1, B), WithinRel(region, 1.0e-12));
}

}  // namespace palace
