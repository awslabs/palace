// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <complex>
#include <vector>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "fem/bilinearform.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/units.hpp"

namespace palace
{
using namespace Catch;

namespace
{

class PaOrderThresholdGuard
{
private:
  int threshold;

public:
  PaOrderThresholdGuard() : threshold(BilinearForm::pa_order_threshold) {}
  ~PaOrderThresholdGuard() { BilinearForm::pa_order_threshold = threshold; }
  PaOrderThresholdGuard(const PaOrderThresholdGuard &) = delete;
  PaOrderThresholdGuard &operator=(const PaOrderThresholdGuard &) = delete;
};

}  // namespace

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

  SECTION("Pole support uses global material indices")
  {
    material2.permittivity_poles.push_back({{-2.0, 5.0}, {1.2, -0.4}});
    MaterialOperator mat_op({material1, material2}, periodic, ProblemType::DRIVEN,
                            palace_mesh);
    CHECK(mat_op.HasPermittivityPoles());
    REQUIRE(mat_op.NumPermittivityPoleMaterials() == 2);
    CHECK(mat_op.HasNonzeroPermittivityPoles(1));
    CHECK(mat_op.HasPermittivityPoleSupport(1));
    CHECK(mat_op.GetPermittivityPoleAttributes(1) == std::vector<int>{2});
    const std::complex<double> s{0.4, 1.1};
    const std::complex<double> pole{-2.0, 5.0}, residue{1.2, -0.4};
    const auto expected =
        residue * s * s / (s - pole) + std::conj(residue) * s * s / (s - std::conj(pole));
    CHECK(std::abs(mat_op.EvaluatePermittivityPoleA2(1, s) - expected) == Approx(0.0));
  }
}

TEST_CASE("MaterialOperator scalar permittivity pole identities",
          "[materialoperator][Serial]")
{
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::TETRAHEDRON));
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh palace_mesh(std::move(par_mesh));
  config::PeriodicBoundaryData periodic;

  config::MaterialData material;
  material.attributes = {1};
  constexpr double eps_inf = 2.5, wp2 = 9.0, gamma = 3.0;
  material.epsilon_r = config::SymmetricMatrixData<3>(eps_inf);
  material.permittivity_poles.push_back({{0.0, 0.0}, {wp2 / gamma, 0.0}});
  material.permittivity_poles.push_back({{-gamma, 0.0}, {-wp2 / gamma, 0.0}});

  MaterialOperator mat_op({material}, periodic, ProblemType::DRIVEN, palace_mesh);
  const std::complex<double> s{1.25, 4.0};
  const std::complex<double> lhs = eps_inf * s * s + mat_op.GetConductivity(1)(0, 0) * s +
                                   mat_op.EvaluatePermittivityPoleA2(0, s);
  const std::complex<double> rhs = eps_inf * s * s + wp2 * s / (s + gamma);
  CHECK(lhs.real() == Approx(rhs.real()));
  CHECK(lhs.imag() == Approx(rhs.imag()));
  CHECK(mat_op.GetConductivity(1)(0, 0) == Approx(wp2 / gamma));
  CHECK(mat_op.HasConductivity());
  CHECK(mat_op.HasPermittivityPoleA2());
  CHECK_THROWS_WITH(mat_op.EvaluatePermittivityPoleA2(0, {-gamma, 0.0}),
                    Catch::Matchers::ContainsSubstring("singular at s = pole"));

  // This distant-pole value is O(1). The forbidden long-division K/C/A2 split forms
  // O(1e16) terms whose cancellation cannot reproduce this result in binary64.
  config::MaterialData distant;
  distant.attributes = {1};
  distant.permittivity_poles.push_back({{-1.0e8, 0.0}, {1.0e8, 0.0}});
  MaterialOperator distant_op({distant}, periodic, ProblemType::DRIVEN, palace_mesh);
  const std::complex<double> distant_expected = -1.0e8 / std::complex<double>(1.0e8, 1.0);
  CHECK(std::abs(distant_op.EvaluatePermittivityPoleA2(0, {0.0, 1.0}) - distant_expected) <
        2.0e-16);

  auto anisotropic = material;
  anisotropic.epsilon_r.s = {eps_inf, 2.0 * eps_inf, eps_inf};
  CHECK_THROWS(MaterialOperator({anisotropic}, periodic, ProblemType::DRIVEN, palace_mesh));
  auto lossy = material;
  lossy.tandelta = config::SymmetricMatrixData<3>(0.1);
  CHECK_THROWS(MaterialOperator({lossy}, periodic, ProblemType::DRIVEN, palace_mesh));
  auto conducting = material;
  conducting.sigma = config::SymmetricMatrixData<3>(1.0);
  MaterialOperator conducting_op({conducting}, periodic, ProblemType::DRIVEN, palace_mesh);
  CHECK(conducting_op.GetConductivity(1)(0, 0) == Approx(1.0 + wp2 / gamma));
}

TEST_CASE("SpaceOperator direct dispersive A2 action", "[materialoperator][Serial]")
{
  auto MakeMesh = []
  {
    auto serial_mesh = std::make_unique<mfem::Mesh>(
        mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::TETRAHEDRON));
    auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(std::move(par_mesh)));
    return mesh;
  };

  config::SolverData solver;
  solver.order = 1;
  solver.linear.mg_max_levels = 1;
  solver.linear.pc_mat_shifted = 0;
  fem::DefaultIntegrationOrder::p_trial = solver.order;
  fem::DefaultIntegrationOrder::q_order_jac = solver.q_order_jac;
  fem::DefaultIntegrationOrder::q_order_extra_pk = solver.q_order_extra;
  fem::DefaultIntegrationOrder::q_order_extra_qk = solver.q_order_extra;
  config::BoundaryData boundaries;
  Units units(1.0, 1.0);

  config::MaterialData dispersive_material;
  dispersive_material.attributes = {1};
  dispersive_material.epsilon_r = config::SymmetricMatrixData<3>(1.0);

  config::MaterialData reference_material;
  reference_material.attributes = {1};
  reference_material.epsilon_r = config::SymmetricMatrixData<3>(1.0);
  config::DomainData reference_domains;
  reference_domains.materials = {reference_material};

  PaOrderThresholdGuard threshold_guard;
  for (int threshold : {1, 2})
  {
    BilinearForm::pa_order_threshold = threshold;  // Partial, then full assembly.
    for (bool distant_pole : {false, true})
    {
      if (distant_pole)
      {
        // This contribution is O(1) near s = i, while the obsolete long-division terms
        // are O(1e16). Exercise the complete operator action, not only scalar evaluation.
        dispersive_material.permittivity_poles = {{{-1.0e8, 0.0}, {1.0e8, 0.0}}};
      }
      else
      {
        // The complete nonzero-pole contribution cancels at s = i but not at s = 2i.
        // This guards structural nonlinear detection against target-frequency sampling.
        dispersive_material.permittivity_poles = {{{-1.0, 0.0}, {4.0, 0.0}},
                                                  {{-2.0, 0.0}, {-20.0, 0.0}},
                                                  {{-3.0, 0.0}, {20.0, 0.0}}};
      }
      config::DomainData dispersive_domains;
      dispersive_domains.materials = {dispersive_material};

      for (bool pc_mat_real : {false, true})
      {
        solver.linear.pc_mat_real = pc_mat_real;
        auto dispersive_mesh = MakeMesh();
        auto reference_mesh = MakeMesh();
        SpaceOperator dispersive(solver, dispersive_domains, boundaries,
                                 ProblemType::DRIVEN, units, dispersive_mesh);
        SpaceOperator reference(solver, reference_domains, boundaries, ProblemType::DRIVEN,
                                units, reference_mesh);

        auto K = dispersive.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
        auto C = dispersive.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
        auto M = dispersive.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
        const std::complex<double> omega{1.3, -0.4};
        const std::complex<double> s = std::complex<double>(0.0, 1.0) * omega;
        auto A2 = dispersive.GetExtraSystemMatrix(omega, Operator::DIAG_ZERO);
        auto K0 = reference.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
        auto B = reference.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
        REQUIRE(K);
        CHECK_FALSE(C);
        REQUIRE(M);
        REQUIRE(A2);
        REQUIRE(K0);
        REQUIRE(B);

        ComplexVector x(K->Width()), y(K->Height()), expected(K->Height()),
            tmp(K->Height());
        x.UseDevice(true);
        y.UseDevice(true);
        expected.UseDevice(true);
        tmp.UseDevice(true);
        linalg::SetRandom(Mpi::World(), x);

        auto CheckAction = [&](const ComplexOperator &op, std::complex<double> coeff,
                               const ComplexOperator *base = nullptr)
        {
          op.Mult(x, y);
          B->Mult(x, expected);
          expected *= coeff;
          if (base)
          {
            base->Mult(x, tmp);
            expected += tmp;
          }
          y.Add(-1.0, expected);
          CHECK(linalg::Norml2(Mpi::World(), y) <
                1.0e-11 * std::max(1.0, linalg::Norml2(Mpi::World(), expected)));
        };

        const auto &mat = dispersive.GetMaterialOp();
        CHECK(mat.HasPermittivityPoleA2());
        if (!distant_pole)
        {
          CHECK(std::abs(mat.EvaluatePermittivityPoleA2(0, {0.0, 1.0})) < 1.0e-14);
          CHECK(std::abs(mat.EvaluatePermittivityPoleA2(0, {0.0, 2.0})) > 1.0e-3);
          auto A2_cancelled =
              dispersive.GetExtraSystemMatrix<ComplexOperator>(1.0, Operator::DIAG_ZERO);
          REQUIRE(A2_cancelled);
          REQUIRE(A2_cancelled->Real());
          REQUIRE(A2_cancelled->Imag());
          A2_cancelled->Mult(x, y);
          CHECK(linalg::Norml2(Mpi::World(), y) < 1.0e-12);
        }
        CheckAction(*K, 0.0, K0.get());  // No pole-derived stiffness contribution.
        CheckAction(*M, 1.0);            // epsilon_inf remains in ordinary M.
        std::complex<double> g = 0.0;
        for (const auto &term : dispersive_material.permittivity_poles)
        {
          g += term.residue * s * s / (s - term.pole);
        }
        CHECK(std::abs(mat.EvaluatePermittivityPoleA2(0, s) - g) <
              1.0e-14 * std::max(1.0, std::abs(g)));
        CheckAction(*A2, g);

        auto system = dispersive.GetSystemMatrix(std::complex<double>{1.0}, s, s * s,
                                                 K.get(), C.get(), M.get(), A2.get());
        K0->Mult(x, expected);
        B->Mult(x, tmp);
        tmp *= s * s + g;
        expected += tmp;
        system->Mult(x, y);
        y.Add(-1.0, expected);
        CHECK(linalg::Norml2(Mpi::World(), y) <
              1.0e-11 * std::max(1.0, linalg::Norml2(Mpi::World(), expected)));

        auto pc = dispersive.GetPreconditionerMatrix<ComplexOperator>(
            std::complex<double>{1.0}, s, s * s, omega);
        const auto *mg =
            dynamic_cast<const BaseMultigridOperator<ComplexOperator> *>(pc.get());
        REQUIRE(mg);
        if (pc_mat_real)
        {
          K0->Mult(x, expected);
          B->Mult(x, tmp);
          tmp *= (s * s).real() + g.real() + g.imag();
          expected += tmp;
        }
        mg->GetFinestOperator().Mult(x, y);
        y.Add(-1.0, expected);
        CHECK(linalg::Norml2(Mpi::World(), y) <
              1.0e-11 * std::max(1.0, linalg::Norml2(Mpi::World(), expected)));
      }
    }
  }
}

TEST_CASE("SpaceOperator zero permittivity poles are ordinary conductivity",
          "[materialoperator][Serial]")
{
  auto MakeMesh = []
  {
    auto serial_mesh = std::make_unique<mfem::Mesh>(
        mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::TETRAHEDRON));
    auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(std::move(par_mesh)));
    return mesh;
  };

  config::SolverData solver;
  solver.order = 1;
  solver.linear.mg_max_levels = 1;
  fem::DefaultIntegrationOrder::p_trial = solver.order;
  fem::DefaultIntegrationOrder::q_order_jac = solver.q_order_jac;
  fem::DefaultIntegrationOrder::q_order_extra_pk = solver.q_order_extra;
  fem::DefaultIntegrationOrder::q_order_extra_qk = solver.q_order_extra;

  config::MaterialData pole_material;
  pole_material.attributes = {1};
  pole_material.epsilon_r = config::SymmetricMatrixData<3>(2.0);
  pole_material.sigma = config::SymmetricMatrixData<3>(1.0);
  pole_material.permittivity_poles.push_back({{0.0, 0.0}, {3.0, 0.0}});
  config::DomainData pole_domains;
  pole_domains.materials = {pole_material};

  auto conductivity_material = pole_material;
  conductivity_material.permittivity_poles.clear();
  conductivity_material.sigma = config::SymmetricMatrixData<3>(4.0);
  config::DomainData conductivity_domains;
  conductivity_domains.materials = {conductivity_material};

  config::BoundaryData boundaries;
  Units units(1.0, 1.0);
  PaOrderThresholdGuard threshold_guard;
  BilinearForm::pa_order_threshold = 1;
  auto pole_mesh = MakeMesh();
  auto conductivity_mesh = MakeMesh();
  SpaceOperator pole(solver, pole_domains, boundaries, ProblemType::DRIVEN, units,
                     pole_mesh);
  SpaceOperator conductivity(solver, conductivity_domains, boundaries, ProblemType::DRIVEN,
                             units, conductivity_mesh);

  const auto &mat = pole.GetMaterialOp();
  CHECK(mat.HasPermittivityPoles());
  CHECK_FALSE(mat.HasNonzeroPermittivityPoles(0));
  CHECK_FALSE(mat.HasPermittivityPoleA2());
  CHECK(mat.HasConductivity());
  CHECK(mat.GetConductivity(1)(0, 0) == Approx(4.0));
  CHECK_FALSE(pole.GetExtraSystemMatrix<ComplexOperator>(1.0, Operator::DIAG_ZERO));

  auto C = pole.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto C_ref = conductivity.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  REQUIRE(C);
  REQUIRE(C_ref);
  ComplexVector x(C->Width()), y(C->Height()), expected(C->Height());
  x.UseDevice(true);
  y.UseDevice(true);
  expected.UseDevice(true);
  linalg::SetRandom(Mpi::World(), x);
  C->Mult(x, y);
  C_ref->Mult(x, expected);
  y.Add(-1.0, expected);
  CHECK(linalg::Norml2(Mpi::World(), y) <
        1.0e-12 * std::max(1.0, linalg::Norml2(Mpi::World(), expected)));
}

TEST_CASE("SpaceOperator dispersive mass action across partitioned material support",
          "[materialoperator][Parallel]")
{
  const MPI_Comm comm = Mpi::World();
  const int size = Mpi::Size(comm);

  auto MakeMesh = [comm, size]
  {
    // Attribute 2 is confined to rank 0's first cell, separated from the rank interface by
    // attribute 1. This reproduces the no-local-support rank condition for the cached B_m.
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
    auto par_mesh =
        std::make_unique<mfem::ParMesh>(comm, *serial_mesh, partitioning.data());
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(std::move(par_mesh)));
    return mesh;
  };

  config::SolverData solver;
  solver.order = 1;
  solver.linear.mg_max_levels = 1;
  solver.linear.pc_mat_real = false;
  solver.linear.pc_mat_shifted = 0;
  fem::DefaultIntegrationOrder::p_trial = solver.order;
  fem::DefaultIntegrationOrder::q_order_jac = solver.q_order_jac;
  fem::DefaultIntegrationOrder::q_order_extra_pk = solver.q_order_extra;
  fem::DefaultIntegrationOrder::q_order_extra_qk = solver.q_order_extra;
  PaOrderThresholdGuard threshold_guard;
  BilinearForm::pa_order_threshold = 1;  // Exercise the partially assembled cache path.

  config::MaterialData material1;
  material1.attributes = {1};
  material1.epsilon_r = config::SymmetricMatrixData<3>(1.0);
  config::MaterialData material2;
  material2.attributes = {2};
  material2.epsilon_r = config::SymmetricMatrixData<3>(1.0);
  material2.permittivity_poles.push_back({{-2.0, 0.0}, {3.0, 0.0}});
  config::DomainData dispersive_domains;
  dispersive_domains.materials = {material1, material2};

  auto reference_material2 = material2;
  reference_material2.permittivity_poles.clear();
  config::DomainData reference_domains;
  reference_domains.materials = {material1, reference_material2};
  auto doubled_material2 = reference_material2;
  doubled_material2.epsilon_r = config::SymmetricMatrixData<3>(2.0);
  config::DomainData doubled_domains;
  doubled_domains.materials = {material1, doubled_material2};

  config::BoundaryData boundaries;
  Units units(1.0, 1.0);
  auto dispersive_mesh = MakeMesh();
  auto reference_mesh = MakeMesh();
  auto doubled_mesh = MakeMesh();
  SpaceOperator dispersive(solver, dispersive_domains, boundaries, ProblemType::DRIVEN,
                           units, dispersive_mesh);
  SpaceOperator reference(solver, reference_domains, boundaries, ProblemType::DRIVEN, units,
                          reference_mesh);
  SpaceOperator doubled(solver, doubled_domains, boundaries, ProblemType::DRIVEN, units,
                        doubled_mesh);
  CHECK((dispersive.GetMaterialOp().GetCeedAttributes(std::vector<int>{2}).Size() > 0) ==
        (Mpi::Rank(comm) == 0));
  CHECK(dispersive.GetMaterialOp().HasPermittivityPoleSupport(1));

  auto K = dispersive.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto C = dispersive.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M = dispersive.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  const std::complex<double> omega{1.3, -0.4};
  auto A2 = dispersive.GetExtraSystemMatrix(omega, Operator::DIAG_ZERO);
  auto K_ref = reference.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M_ref = reference.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M_doubled = doubled.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  REQUIRE(K);
  CHECK_FALSE(C);
  REQUIRE(M);
  REQUIRE(A2);
  REQUIRE(K_ref);
  REQUIRE(M_ref);
  REQUIRE(M_doubled);

  ComplexVector x(K->Width()), y(K->Height()), expected(K->Height()), tmp(K->Height()),
      B2x(K->Height());
  x.UseDevice(true);
  y.UseDevice(true);
  expected.UseDevice(true);
  tmp.UseDevice(true);
  B2x.UseDevice(true);
  linalg::SetRandom(comm, x);

  auto ApplyB2 = [&]
  {
    M_doubled->Mult(x, B2x);
    M_ref->Mult(x, tmp);
    B2x.Add(-1.0, tmp);
  };
  auto CheckAction = [&](const ComplexOperator &op, std::complex<double> coeff,
                         const ComplexOperator *base = nullptr)
  {
    ApplyB2();
    expected = B2x;
    expected *= coeff;
    if (base)
    {
      base->Mult(x, tmp);
      expected += tmp;
    }
    op.Mult(x, y);
    y.Add(-1.0, expected);
    CHECK(linalg::Norml2(comm, y) <
          1.0e-11 * std::max(1.0, linalg::Norml2(comm, expected)));
  };

  const auto &mat = dispersive.GetMaterialOp();
  const std::complex<double> s = std::complex<double>(0.0, 1.0) * omega;
  CheckAction(*K, 0.0, K_ref.get());
  M_ref->Mult(x, expected);
  M->Mult(x, y);
  y.Add(-1.0, expected);
  CHECK(linalg::Norml2(comm, y) < 1.0e-11 * std::max(1.0, linalg::Norml2(comm, expected)));
  CheckAction(*A2, mat.EvaluatePermittivityPoleA2(1, s));

  auto system = dispersive.GetSystemMatrix(std::complex<double>{1.0}, s, s * s, K.get(),
                                           C.get(), M.get(), A2.get());
  auto pc = dispersive.GetPreconditionerMatrix<ComplexOperator>(std::complex<double>{1.0},
                                                                s, s * s, omega);
  const auto *mg = dynamic_cast<const BaseMultigridOperator<ComplexOperator> *>(pc.get());
  REQUIRE(mg);
  system->Mult(x, expected);
  mg->GetFinestOperator().Mult(x, y);
  y.Add(-1.0, expected);
  CHECK(linalg::Norml2(comm, y) < 1.0e-11 * std::max(1.0, linalg::Norml2(comm, expected)));
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

}  // namespace palace
