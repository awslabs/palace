// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <memory>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include "fem/coefficient.hpp"
#include "fem/domain_point_field_evaluator.hpp"
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/labels.hpp"

namespace palace
{

namespace
{

void CheckField(const mfem::ParGridFunction &value, const mfem::ParGridFunction &reference)
{
  const double *v = value.HostRead();
  const double *r = reference.HostRead();
  double max_diff = 0.0, max_ref = 0.0;
  for (int i = 0; i < reference.Size(); i++)
  {
    max_diff = std::max(max_diff, std::abs(v[i] - r[i]));
    max_ref = std::max(max_ref, std::abs(r[i]));
  }
  CAPTURE(max_diff, max_ref);
  CHECK(max_diff <= 1.0e-11 * std::max(max_ref, 1.0));
}

}  // namespace

TEST_CASE("Domain point field evaluator 3D", "[postoperator][Serial][Parallel][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const auto elem_type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  const auto order = GENERATE(1, 2);
  const auto complex = GENERATE(false, true);
  CAPTURE(elem_type, order, complex);
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  auto smesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(2, 2, 2, elem_type, 1.2, 0.9, 0.7));
  smesh->EnsureNodes();
  REQUIRE(Mpi::Size(comm) <= smesh->GetNE());
  auto pmesh = std::make_unique<mfem::ParMesh>(comm, *smesh);
  Mesh mesh(std::move(pmesh));

  config::MaterialData material;
  material.attributes = {1};
  material.epsilon_r.s[0] = 2.0;
  material.epsilon_r.s[1] = 3.0;
  material.epsilon_r.s[2] = 4.0;
  material.mu_r.s[0] = 1.0;
  material.mu_r.s[1] = 1.5;
  material.mu_r.s[2] = 2.0;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({material}, periodic, ProblemType::DRIVEN, mesh);

  mfem::ND_FECollection nd_fec(order, 3);
  mfem::RT_FECollection rt_fec(order - 1, 3);
  FiniteElementSpace nd_fespace(mesh, &nd_fec), rt_fespace(mesh, &rt_fec);
  GridFunction E(nd_fespace, complex), B(rt_fespace, complex);

  mfem::VectorFunctionCoefficient e_real(3,
                                         [](const mfem::Vector &x, mfem::Vector &v)
                                         {
                                           v(0) = std::sin(x(1)) + x(2) * x(2);
                                           v(1) = std::cos(x(2)) + x(0);
                                           v(2) = x(0) * x(1) + 1.0;
                                         });
  mfem::VectorFunctionCoefficient b_real(3,
                                         [](const mfem::Vector &x, mfem::Vector &v)
                                         {
                                           v(0) = x(1) - 0.3 * x(2);
                                           v(1) = std::sin(x(2)) + 0.5;
                                           v(2) = std::cos(x(0)) - x(1) * x(2);
                                         });
  E.Real().ProjectCoefficient(e_real);
  B.Real().ProjectCoefficient(b_real);
  if (complex)
  {
    mfem::VectorFunctionCoefficient e_imag(3,
                                           [](const mfem::Vector &x, mfem::Vector &v)
                                           {
                                             v(0) = x(1) * x(2) - 0.5;
                                             v(1) = std::sin(x(0)) - x(2);
                                             v(2) = std::cos(x(1)) + x(0) * x(0);
                                           });
    mfem::VectorFunctionCoefficient b_imag(3,
                                           [](const mfem::Vector &x, mfem::Vector &v)
                                           {
                                             v(0) = std::cos(x(2)) - 0.2;
                                             v(1) = x(0) * x(2) + 0.1;
                                             v(2) = std::sin(x(1)) - x(0);
                                           });
    E.Imag().ProjectCoefficient(e_imag);
    B.Imag().ProjectCoefficient(b_imag);
  }

  mfem::L2_FECollection output_fec(order, 3);
  mfem::ParFiniteElementSpace scalar_fespace(&mesh.Get(), &output_fec);
  mfem::ParFiniteElementSpace vector_fespace(&mesh.Get(), &output_fec, 3);
  const double scaling = 2.5;

  SECTION("Electric energy")
  {
    DomainPointFieldEvaluator evaluator(DomainPointFieldEvaluator::Kind::ENERGY_E, mesh,
                                        mat_op, E.ParFESpace(), nullptr, scalar_fespace,
                                        scaling);
    REQUIRE(evaluator.IsValid());
    mfem::ParGridFunction value(&scalar_fespace), reference(&scalar_fespace);
    evaluator.Eval(&E, nullptr, value);
    EnergyDensityCoefficient<EnergyDensityType::ELECTRIC> legacy(E, mat_op, scaling);
    reference.ProjectCoefficient(legacy);
    CheckField(value, reference);
  }

  SECTION("Magnetic energy")
  {
    DomainPointFieldEvaluator evaluator(DomainPointFieldEvaluator::Kind::ENERGY_M, mesh,
                                        mat_op, nullptr, B.ParFESpace(), scalar_fespace,
                                        scaling);
    REQUIRE(evaluator.IsValid());
    mfem::ParGridFunction value(&scalar_fespace), reference(&scalar_fespace);
    evaluator.Eval(nullptr, &B, value);
    EnergyDensityCoefficient<EnergyDensityType::MAGNETIC> legacy(B, mat_op, scaling);
    reference.ProjectCoefficient(legacy);
    CheckField(value, reference);
  }

  SECTION("Poynting vector")
  {
    DomainPointFieldEvaluator evaluator(DomainPointFieldEvaluator::Kind::POYNTING, mesh,
                                        mat_op, E.ParFESpace(), B.ParFESpace(),
                                        vector_fespace, scaling);
    REQUIRE(evaluator.IsValid());
    mfem::ParGridFunction value(&vector_fespace), reference(&vector_fespace);
    evaluator.Eval(&E, &B, value);
    PoyntingVectorCoefficient legacy(E, B, mat_op, scaling);
    reference.ProjectCoefficient(legacy);
    CheckField(value, reference);
  }
}

TEST_CASE("Domain point field evaluator 2D", "[postoperator][Serial][Parallel][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const auto elem_type = GENERATE(mfem::Element::TRIANGLE, mfem::Element::QUADRILATERAL);
  const auto order = GENERATE(1, 2);
  const auto complex = GENERATE(false, true);
  CAPTURE(elem_type, order, complex);
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  auto smesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(3, 2, elem_type, false, 1.2, 0.7));
  smesh->EnsureNodes();
  REQUIRE(Mpi::Size(comm) <= smesh->GetNE());
  auto pmesh = std::make_unique<mfem::ParMesh>(comm, *smesh);
  Mesh mesh(std::move(pmesh));

  config::MaterialData material;
  material.attributes = {1};
  material.epsilon_r.s[0] = 2.0;
  material.epsilon_r.s[1] = 3.0;
  material.epsilon_r.s[2] = 4.0;
  material.mu_r.s[0] = 1.0;
  material.mu_r.s[1] = 1.5;
  material.mu_r.s[2] = 2.0;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({material}, periodic, ProblemType::DRIVEN, mesh);

  mfem::ND_FECollection nd_fec(order, 2);
  mfem::L2_FECollection l2_fec(order - 1, 2, mfem::BasisType::GaussLegendre,
                               mfem::FiniteElement::INTEGRAL);
  FiniteElementSpace nd_fespace(mesh, &nd_fec), l2_fespace(mesh, &l2_fec);
  GridFunction E(nd_fespace, complex), B(l2_fespace, complex);

  mfem::VectorFunctionCoefficient e_real(2,
                                         [](const mfem::Vector &x, mfem::Vector &v)
                                         {
                                           v(0) = std::sin(x(1)) + 0.2 * x(0);
                                           v(1) = std::cos(x(0)) + x(1) * x(1);
                                         });
  mfem::FunctionCoefficient b_real([](const mfem::Vector &x)
                                   { return 0.3 + x(0) - 0.7 * x(1); });
  E.Real().ProjectCoefficient(e_real);
  B.Real().ProjectCoefficient(b_real);
  if (complex)
  {
    mfem::VectorFunctionCoefficient e_imag(2,
                                           [](const mfem::Vector &x, mfem::Vector &v)
                                           {
                                             v(0) = x(0) * x(1) - 0.1;
                                             v(1) = std::sin(x(0) + x(1));
                                           });
    mfem::FunctionCoefficient b_imag([](const mfem::Vector &x)
                                     { return std::cos(x(0)) + 0.4 * x(1); });
    E.Imag().ProjectCoefficient(e_imag);
    B.Imag().ProjectCoefficient(b_imag);
  }

  mfem::L2_FECollection output_fec(order, 2);
  mfem::ParFiniteElementSpace scalar_fespace(&mesh.Get(), &output_fec);
  mfem::ParFiniteElementSpace vector_fespace(&mesh.Get(), &output_fec, 2);
  const double scaling = 1.7;

  SECTION("Electric energy")
  {
    DomainPointFieldEvaluator evaluator(DomainPointFieldEvaluator::Kind::ENERGY_E, mesh,
                                        mat_op, E.ParFESpace(), nullptr, scalar_fespace,
                                        scaling);
    REQUIRE(evaluator.IsValid());
    mfem::ParGridFunction value(&scalar_fespace), reference(&scalar_fespace);
    evaluator.Eval(&E, nullptr, value);
    EnergyDensityCoefficient<EnergyDensityType::ELECTRIC> legacy(E, mat_op, scaling);
    reference.ProjectCoefficient(legacy);
    CheckField(value, reference);
  }

  SECTION("Magnetic energy")
  {
    DomainPointFieldEvaluator evaluator(DomainPointFieldEvaluator::Kind::ENERGY_M, mesh,
                                        mat_op, nullptr, B.ParFESpace(), scalar_fespace,
                                        scaling);
    REQUIRE(evaluator.IsValid());
    mfem::ParGridFunction value(&scalar_fespace), reference(&scalar_fespace);
    evaluator.Eval(nullptr, &B, value);
    EnergyDensityCoefficient<EnergyDensityType::MAGNETIC> legacy(B, mat_op, scaling);
    reference.ProjectCoefficient(legacy);
    CheckField(value, reference);
  }

  SECTION("Poynting vector")
  {
    DomainPointFieldEvaluator evaluator(DomainPointFieldEvaluator::Kind::POYNTING, mesh,
                                        mat_op, E.ParFESpace(), B.ParFESpace(),
                                        vector_fespace, scaling);
    REQUIRE(evaluator.IsValid());
    mfem::ParGridFunction value(&vector_fespace), reference(&vector_fespace);
    evaluator.Eval(&E, &B, value);
    PoyntingVectorCoefficient legacy(E, B, mat_op, scaling);
    reference.ProjectCoefficient(legacy);
    CheckField(value, reference);
  }
}

}  // namespace palace
