// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <memory>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "fem/bilinearform.hpp"
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "linalg/vector.hpp"
#include "models/domainpostoperator.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/constants.hpp"
#include "utils/units.hpp"

namespace palace
{
using namespace Catch::Matchers;

namespace
{

// Rectangle [a, b] x [0, h] of the (r, z) half-plane, meshed with triangles.
std::unique_ptr<mfem::ParMesh> MakeAnnulusSection(MPI_Comm comm, double a, double b,
                                                  double h, int nr, int nz)
{
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(nr, nz, mfem::Element::TRIANGLE, true, b - a, h));
  for (int i = 0; i < serial_mesh->GetNV(); i++)
  {
    serial_mesh->GetVertex(i)[0] += a;
  }
  return std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
}

double UniformFieldEnergy(Mesh &palace_mesh, const Units &units, int order)
{
  config::MaterialData material;
  material.attributes = {1};
  std::vector<config::MaterialData> materials = {material};
  config::PeriodicBoundaryData periodic;
  config::DomainPostData postpro;
  MaterialOperator mat_op(materials, periodic, ProblemType::ELECTROSTATIC, palace_mesh);
  mfem::H1_FECollection h1_fec(order, palace_mesh.Dimension());
  FiniteElementSpace h1_fespace(palace_mesh, &h1_fec);
  DomainPostOperator dom_post_op(postpro, mat_op, h1_fespace);

  // V = -E0 z gives the uniform field E = E0 z_hat (SI units: V/m).
  const double E0 = units.Nondimensionalize<Units::ValueType::FIELD_E>(1.0e6);
  auto V_linear = [&](const mfem::Vector &x) { return -E0 * x(1); };
  GridFunction V(h1_fespace, false);
  mfem::FunctionCoefficient V_coeff(V_linear);
  V.Real().ProjectCoefficient(V_coeff);
  return units.Dimensionalize<Units::ValueType::ENERGY>(
      dom_post_op.GetElectricFieldEnergy(V));
}

}  // namespace

TEST_CASE("Axisymmetric - libCEED geometry factors carry 2 pi r", "[axisymmetric][Serial]")
{
  // Uniform axial field in the annular section [a, b] x [0, h]: the Cartesian (translation
  // invariant, depth Lc) energy is 1/2 eps0 E0^2 (b - a) h Lc, the axisymmetric one is
  // 1/2 eps0 E0^2 pi (b^2 - a^2) h. Both are exact for any quadrature of the linear V.
  constexpr int order = 2;
  fem::DefaultIntegrationOrder::p_trial = order;
  constexpr double a = 0.7, b = 2.3, h = 1.9;
  Units units(1.0e-6, 3.0e-6);  // L0 = 1 um, Lc = 3 um
  MPI_Comm comm = Mpi::World();

  Mesh palace_mesh(MakeAnnulusSection(comm, a, b, h, 6, 5));
  REQUIRE_FALSE(palace_mesh.IsAxisymmetric());
  const double cartesian = UniformFieldEnergy(palace_mesh, units, order);

  palace_mesh.SetAxisymmetric(true);
  REQUIRE(palace_mesh.IsAxisymmetric());
  const double axisymmetric = UniformFieldEnergy(palace_mesh, units, order);

  // The test mesh coordinates are already nondimensional (units of Lc), so one mesh unit
  // is Lc = 3 um and the Cartesian depth is Lc.
  const double E0 = 1.0e6;
  const double L = units.Dimensionalize<Units::ValueType::LENGTH>(1.0);
  CHECK_THAT(L, WithinRel(3.0e-6, 1.0e-12));
  const double expected_cartesian =
      0.5 * electromagnetics::epsilon0_ * E0 * E0 * ((b - a) * L) * (h * L) * L;
  const double expected_axisymmetric = 0.5 * electromagnetics::epsilon0_ * E0 * E0 * M_PI *
                                       ((b * b - a * a) * L * L) * (h * L);
  CHECK_THAT(cartesian, WithinRel(expected_cartesian, 1.0e-10));
  CHECK_THAT(axisymmetric, WithinRel(expected_axisymmetric, 1.0e-10));
  CHECK_THAT(axisymmetric / cartesian,
             WithinRel(M_PI * (b + a), 1.0e-10));  // 2 pi r_mean / depth, depth = 1

  // Switching the flag back restores the Cartesian measure exactly.
  palace_mesh.SetAxisymmetric(false);
  CHECK(UniformFieldEnergy(palace_mesh, units, order) == cartesian);
}

TEST_CASE("Axisymmetric - boundary geometry factors carry 2 pi r", "[axisymmetric][Serial]")
{
  // Boundary mass matrix on the whole boundary of [a, b] x [0, h] applied to the constant
  // function: sum = perimeter (Cartesian) or the area of revolution
  // 2 pi [ (b^2 - a^2) + (a + b) h ] (axisymmetric).
  constexpr int order = 1;
  fem::DefaultIntegrationOrder::p_trial = order;
  constexpr double a = 0.5, b = 1.5, h = 2.0;
  MPI_Comm comm = Mpi::World();
  Mesh palace_mesh(MakeAnnulusSection(comm, a, b, h, 4, 8));

  auto BoundaryMeasure = [&]()
  {
    mfem::H1_FECollection h1_fec(order, palace_mesh.Dimension());
    FiniteElementSpace h1_fespace(palace_mesh, &h1_fec);
    std::vector<int> bdr_attributes;
    for (int attr = 1; attr <= palace_mesh.Get().bdr_attributes.Max(); attr++)
    {
      bdr_attributes.push_back(attr);
    }
    MaterialPropertyCoefficient one_func(palace_mesh.MaxCeedBdrAttribute());
    one_func.AddMaterialProperty(palace_mesh.GetCeedBdrAttributes(bdr_attributes), 1.0);
    BilinearForm m(h1_fespace);
    m.AddBoundaryIntegrator<MassIntegrator>(one_func);
    auto M = m.FullAssemble(false);
    Vector x(M->Height()), y(M->Height());
    x = 1.0;
    M->Mult(x, y);
    double sum = linalg::LocalSum(y);
    Mpi::GlobalSum(1, &sum, comm);
    return sum;
  };

  const double cartesian = BoundaryMeasure();
  CHECK_THAT(cartesian, WithinRel(2.0 * (b - a) + 2.0 * h, 1.0e-10));
  palace_mesh.SetAxisymmetric(true);
  const double axisymmetric = BoundaryMeasure();
  CHECK_THAT(axisymmetric,
             WithinRel(2.0 * M_PI * ((b * b - a * a) + (a + b) * h), 1.0e-10));
}

}  // namespace palace
