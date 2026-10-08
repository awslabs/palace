// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <complex>
#include <limits>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"
#include "fem/libceed/operator.hpp"
#include "fem/mesh.hpp"
#include "linalg/hypre.hpp"
#include "linalg/rap.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/labels.hpp"
#include "utils/units.hpp"

namespace palace
{
namespace
{

struct IntegrationSettingsGuard
{
  int pa_order_threshold = BilinearForm::pa_order_threshold;
  int p_trial = fem::DefaultIntegrationOrder::p_trial;
  bool q_order_jac = fem::DefaultIntegrationOrder::q_order_jac;
  int q_order_extra_pk = fem::DefaultIntegrationOrder::q_order_extra_pk;
  int q_order_extra_qk = fem::DefaultIntegrationOrder::q_order_extra_qk;

  ~IntegrationSettingsGuard()
  {
    BilinearForm::pa_order_threshold = pa_order_threshold;
    fem::DefaultIntegrationOrder::p_trial = p_trial;
    fem::DefaultIntegrationOrder::q_order_jac = q_order_jac;
    fem::DefaultIntegrationOrder::q_order_extra_pk = q_order_extra_pk;
    fem::DefaultIntegrationOrder::q_order_extra_qk = q_order_extra_qk;
  }
};

struct ImaginaryHierarchyData
{
  std::vector<HYPRE_Int> coarse_row_offsets;
  std::vector<HYPRE_Int> coarse_columns;
  std::vector<double> coarse_values;
  CeedInt fine_suboperators = 0;
};

ImaginaryHierarchyData InspectImaginaryHierarchy(const ComplexOperator &op)
{
  const auto *mg_op = dynamic_cast<const ComplexMultigridOperator *>(&op);
  REQUIRE(mg_op);
  REQUIRE(mg_op->GetNumLevels() == 2);

  const auto *coarse_op =
      dynamic_cast<const ComplexParOperator *>(&mg_op->GetOperatorAtLevel(0));
  REQUIRE(coarse_op);
  const auto *coarse_imag =
      dynamic_cast<const hypre::HypreCSRMatrix *>(coarse_op->LocalOperator().Imag());
  REQUIRE(coarse_imag);
  hypre_CSRMatrixMigrate(*coarse_imag, HYPRE_MEMORY_HOST);

  ImaginaryHierarchyData data;
  data.coarse_row_offsets.assign(coarse_imag->GetI(),
                                 coarse_imag->GetI() + coarse_imag->Height() + 1);
  data.coarse_columns.assign(coarse_imag->GetJ(), coarse_imag->GetJ() + coarse_imag->NNZ());
  data.coarse_values.assign(coarse_imag->GetData(),
                            coarse_imag->GetData() + coarse_imag->NNZ());

  const auto *fine_op =
      dynamic_cast<const ComplexParOperator *>(&mg_op->GetOperatorAtLevel(1));
  REQUIRE(fine_op);
  const auto *fine_imag =
      dynamic_cast<const ceed::Operator *>(fine_op->LocalOperator().Imag());
  REQUIRE(fine_imag);
  for (std::size_t i = 0; i < fine_imag->Size(); i++)
  {
    CeedInt num_suboperators = 0;
    REQUIRE(CeedOperatorCompositeGetNumSub((*fine_imag)[i], &num_suboperators) == 0);
    data.fine_suboperators += num_suboperators;
  }
  return data;
}

}  // namespace

TEST_CASE("SpaceOperator retains coarse support while omitting fine exact zeros",
          "[spaceoperator][Serial][Parallel]")
{
  using namespace std::complex_literals;

  MPI_Comm comm = Mpi::World();
  mfem::Mesh serial_mesh = mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::TETRAHEDRON);
  while (serial_mesh.GetNE() < Mpi::Size(comm))
  {
    serial_mesh.UniformRefinement();
  }
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(comm, serial_mesh));

  IntegrationSettingsGuard settings_guard;
  config::SolverData solver;
  solver.order = 2;
  solver.pa_order_threshold = 2;
  solver.linear.mg_max_levels = 2;
  solver.linear.mg_coarsening = MultigridCoarsening::LINEAR;
  solver.linear.pc_mat_real = false;
  solver.linear.pc_mat_shifted = 0;
  BilinearForm::pa_order_threshold = solver.pa_order_threshold;
  fem::DefaultIntegrationOrder::p_trial = solver.order;
  fem::DefaultIntegrationOrder::q_order_jac = solver.q_order_jac;
  fem::DefaultIntegrationOrder::q_order_extra_pk = solver.q_order_extra;
  fem::DefaultIntegrationOrder::q_order_extra_qk = solver.q_order_extra;

  config::MaterialData material;
  material.attributes = {1};
  config::DomainData domains;
  domains.attributes = {1};
  domains.materials = {material};
  config::BoundaryData boundaries;
  Units units(1.0, 1.0);
  SpaceOperator space_op(solver, domains, boundaries, ProblemType::EIGENMODE, units, mesh);

  constexpr double omega = 2.0;
  const std::complex<double> lambda_zero = 1i * omega;
  const double epsilon = std::numeric_limits<double>::epsilon();
  const std::complex<double> lambda_tiny = epsilon + 1i * omega;
  auto Assemble = [&space_op](std::complex<double> lambda)
  {
    return space_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, lambda,
                                                             lambda * lambda, lambda / 1i);
  };

  auto zero_pc = Assemble(lambda_zero);
  auto zero_data = InspectImaginaryHierarchy(*zero_pc);
  auto tiny_pc = Assemble(lambda_tiny);
  auto tiny_data = InspectImaginaryHierarchy(*tiny_pc);

  CHECK(zero_data.coarse_row_offsets == tiny_data.coarse_row_offsets);
  CHECK(zero_data.coarse_columns == tiny_data.coarse_columns);
  CHECK(std::ranges::all_of(zero_data.coarse_values,
                            [](double value) { return value == 0.0; }));
  CHECK(std::ranges::any_of(tiny_data.coarse_values,
                            [](double value) { return value != 0.0; }));
  CHECK(zero_data.fine_suboperators == 0);
  CHECK(tiny_data.fine_suboperators > 0);
}

TEST_CASE("CPU preconditioner quadrature data preserves the multigrid operators",
          "[spaceoperator][Serial][Parallel]")
{
  const auto element = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  const auto comm = Mpi::World();
  auto serial_mesh = mfem::Mesh::MakeCartesian3D(2, 2, 2, element);
  serial_mesh.SetCurvature(2);
  serial_mesh.Transform(
      [](const mfem::Vector &x, mfem::Vector &y)
      {
        y = x;
        y(0) += 0.05 * x(1) * x(2);
        y(1) += 0.03 * x(0) * x(0);
      });
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(comm, serial_mesh));

  IntegrationSettingsGuard settings_guard;
  config::SolverData solver;
  solver.order = 3;
  solver.pa_order_threshold = 2;
  solver.linear.mg_max_levels = 3;
  solver.linear.mg_coarsening = MultigridCoarsening::LINEAR;
  solver.linear.pc_mat_shifted = 0;
  BilinearForm::pa_order_threshold = solver.pa_order_threshold;
  fem::DefaultIntegrationOrder::p_trial = solver.order;
  fem::DefaultIntegrationOrder::q_order_jac = solver.q_order_jac;
  fem::DefaultIntegrationOrder::q_order_extra_pk = solver.q_order_extra;
  fem::DefaultIntegrationOrder::q_order_extra_qk = solver.q_order_extra;

  config::MaterialData material;
  material.attributes = {1};
  material.mu_r.s = {1.0, 2.0, 3.0};
  material.epsilon_r.s = {2.0, 4.0, 5.0};
  config::DomainData domains;
  domains.attributes = {1};
  domains.materials = {material};
  config::BoundaryData boundaries;
  Units units(1.0, 1.0);
  SpaceOperator space_op(solver, domains, boundaries, ProblemType::EIGENMODE, units, mesh);
  auto pc = space_op.GetPreconditionerMatrix<Operator>(1.5, 0.0, 2.0, 0.0);
  const auto *mg = dynamic_cast<const MultigridOperator *>(pc.get());
  REQUIRE(mg);
  REQUIRE(mg->GetNumLevels() == 3);

  const auto &mat = space_op.GetMaterialOp();
  MaterialPropertyCoefficient curl(mat.GetAttributeToMaterial(),
                                   mat.GetCurlCurlInvPermeability(), 1.5);
  MaterialPropertyCoefficient mass(mat.GetAttributeToMaterial(), mat.GetPermittivityAbs(),
                                   2.0);
  for (bool auxiliary : {false, true})
  {
    const auto &spaces = auxiliary ? space_op.GetH1Spaces() : space_op.GetNDSpaces();
    for (std::size_t l = 0; l < mg->GetNumLevels(); l++)
    {
      const auto &space = spaces.GetFESpaceAtLevel(l);
      // Independently assemble each level without cached quadrature data. This also
      // checks reuse of the cached tensors between polynomial levels.
      BilinearForm form(space);
      if (auxiliary)
      {
        form.AddDomainIntegrator<DiffusionIntegrator>(mass);
      }
      else
      {
        form.AddDomainIntegrator<CurlCurlMassIntegrator>(curl, mass);
      }
      ParOperator reference(form.Assemble(false), space);
      const auto &actual =
          auxiliary ? mg->GetAuxiliaryOperatorAtLevel(l) : mg->GetOperatorAtLevel(l);
      Vector x(space.GetTrueVSize()), expected(x.Size()), result(x.Size());
      for (int seed : {42, 314159})
      {
        x.Randomize(seed + Mpi::Rank(comm));
        reference.Mult(x, expected);
        actual.Mult(x, result);
        result.Add(-1.0, expected);
        CHECK_THAT(linalg::Norml2(comm, result) / linalg::Norml2(comm, expected),
                   Catch::Matchers::WithinAbs(0.0, 2.0e-13));
      }
    }
  }
}

TEST_CASE("SpaceOperator frequency-dependent PML at complex frequency",
          "[spaceoperator][pml][Serial][Parallel]")
{
  using namespace std::complex_literals;

  // Unit cube with a PML layer z ∈ [0.75, 1] (attribute 2).
  MPI_Comm comm = Mpi::World();
  mfem::Mesh serial_mesh = mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::HEXAHEDRON);
  for (int i = 0; i < serial_mesh.GetNE(); i++)
  {
    mfem::Vector center;
    serial_mesh.GetElementCenter(i, center);
    serial_mesh.SetAttribute(i, (center(2) > 0.75) ? 2 : 1);
  }
  serial_mesh.SetAttributes();
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(comm, serial_mesh));

  IntegrationSettingsGuard settings_guard;
  config::SolverData solver;
  solver.order = 2;
  solver.linear.mg_max_levels = 1;
  solver.linear.pc_mat_real = false;
  solver.linear.pc_mat_shifted = 0;
  fem::DefaultIntegrationOrder::p_trial = solver.order;

  constexpr double omega0 = 2.0;
  auto MakeSpaceOperator = [&](bool frequency_dependent)
  {
    config::MaterialData vacuum, pml;
    vacuum.attributes = {1};
    pml.attributes = {2};
    pml.epsilon_r.s = {2.0, 2.0, 2.0};
    pml.pml = config::PMLData();
    pml.pml->autodetect_geometry = false;
    pml.pml->direction_signs = {0, 0, 0, 0, 0, 1};
    pml.pml->thickness = {0.0, 0.0, 0.0, 0.0, 0.0, 0.25};
    pml.pml->kappa_max = {1.0, 1.0, 1.5};
    pml.pml->alpha_max = {0.0, 0.0, 0.2};
    pml.pml->frequency_dependent = frequency_dependent;
    pml.pml->reference_frequency = omega0;
    config::DomainData domains;
    domains.attributes = {1, 2};
    domains.materials = {vacuum, pml};
    config::BoundaryData boundaries;
    Units units(1.0, 1.0);
    return std::make_unique<SpaceOperator>(solver, domains, boundaries,
                                           ProblemType::EIGENMODE, units, mesh);
  };
  auto fd_op = MakeSpaceOperator(true);
  const ComplexOperator *C0 = nullptr;
  REQUIRE(fd_op->GetMaterialOp().HasPML());
  REQUIRE(fd_op->GetMaterialOp().HasFrequencyDependentPML());

  ComplexVector x(fd_op->GetNDSpace().GetTrueVSize()), y1(x.Size()), y2(x.Size());
  x.UseDevice(true);
  y1.UseDevice(true);
  y2.UseDevice(true);
  x.Real().Randomize(42 + Mpi::Rank(comm));
  x.Imag().Randomize(314159 + Mpi::Rank(comm));
  auto RelErr = [&](ComplexVector &a, const ComplexVector &b)
  {
    a.Add(-1.0, b);
    return linalg::Norml2(comm, a) / linalg::Norml2(comm, b);
  };
  auto A2 = [&](std::complex<double> omega)
  {
    auto op = fd_op->GetExtraSystemMatrix(omega, Operator::DIAG_ZERO);
    REQUIRE(op);
    return op;
  };

  SECTION("Complex frequency overload matches the real frequency one for real ω")
  {
    auto A2r = fd_op->GetExtraSystemMatrix<ComplexOperator>(omega0, Operator::DIAG_ZERO);
    REQUIRE(A2r);
    A2r->Mult(x, y1);
    A2(omega0)->Mult(x, y2);
    CHECK_THAT(RelErr(y2, y1), Catch::Matchers::WithinAbs(0.0, 1.0e-14));
  }

  SECTION("System matrix at real ω matches the static PML with ω₀ = ω")
  {
    auto static_op = MakeSpaceOperator(false);
    REQUIRE(!static_op->GetMaterialOp().HasFrequencyDependentPML());
    auto Ks = static_op->GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto Ms = static_op->GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto As = static_op->GetSystemMatrix(1.0 + 0.0i, 0.0 + 0.0i, -omega0 * omega0 + 0.0i,
                                         Ks.get(), C0, Ms.get());
    auto Kf = fd_op->GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto Mf = fd_op->GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto A2f = A2(omega0);
    auto Af = fd_op->GetSystemMatrix(1.0 + 0.0i, 0.0 + 0.0i, -omega0 * omega0 + 0.0i,
                                     Kf.get(), C0, Mf.get(), A2f.get());
    As->Mult(x, y1);
    Af->Mult(x, y2);
    CHECK_THAT(RelErr(y2, y1), Catch::Matchers::WithinAbs(0.0, 1.0e-13));
  }

  SECTION("Frozen-stretch PML matrices are the static PML matrices")
  {
    // The frozen-stretch matrices at ω₀ reproduce A2(ω₀) and the static PML terms of the
    // stiffness and mass matrices for the reference frequency ω₀.
    auto P = fd_op->GetFrequencyDependentPMLMatrices(omega0, Operator::DIAG_ZERO);
    REQUIRE(P[0]);
    REQUIRE(!P[1]);
    REQUIRE(P[2]);
    ComplexVector t(x.Size());
    t.UseDevice(true);
    P[0]->Mult(x, y1);
    P[2]->Mult(x, t);
    y1.Add(-omega0 * omega0, t);
    A2(omega0)->Mult(x, y2);
    CHECK_THAT(RelErr(y1, y2), Catch::Matchers::WithinAbs(0.0, 1.0e-14));

    auto static_op = MakeSpaceOperator(false);
    auto Ks = static_op->GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto Kf = fd_op->GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    Ks->Mult(x, y1);
    Kf->Mult(x, y2);
    P[0]->AddMult(x, y2);
    CHECK_THAT(RelErr(y2, y1), Catch::Matchers::WithinAbs(0.0, 1.0e-14));
    auto Ms = static_op->GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto Mf = fd_op->GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    Ms->Mult(x, y1);
    Mf->Mult(x, y2);
    P[2]->AddMult(x, y2);
    CHECK_THAT(RelErr(y2, y1), Catch::Matchers::WithinAbs(0.0, 1.0e-14));

    // Without the PML terms, A2 is empty.
    CHECK(!fd_op->GetExtraSystemMatrix(std::complex<double>(omega0, 0.3),
                                       Operator::DIAG_ZERO, false));
  }

  SECTION("A2(ω) is the analytic continuation in ω")
  {
    // Cauchy-Riemann equations for the operator action, with central differences.
    const std::complex<double> omega = {omega0, 0.3};
    const double h = 1.0e-5;
    ComplexVector t(x.Size());
    t.UseDevice(true);
    A2(omega + h)->Mult(x, y1);
    A2(omega - h)->Mult(x, t);
    y1.Add(-1.0, t);  // ∂A2/∂Re{ω} x (2 h)
    A2(omega + 1i * h)->Mult(x, y2);
    A2(omega - 1i * h)->Mult(x, t);
    y2.Add(-1.0, t);
    y2 *= -1i;  // -i ∂A2/∂Im{ω} x (2 h) = ∂A2/∂Re{ω} x (2 h)
    CHECK_THAT(RelErr(y2, y1), Catch::Matchers::WithinAbs(0.0, 1.0e-8));

    // The stretch is evaluated at the complex frequency, not its real part.
    A2(omega)->Mult(x, y1);
    A2(omega.real())->Mult(x, y2);
    CHECK(RelErr(y2, y1) > 1.0e-3);
  }

  SECTION("Preconditioner matches the system matrix at complex ω")
  {
    const std::complex<double> omega = {omega0, 0.3}, lambda = 1i * omega;
    auto K = fd_op->GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto M = fd_op->GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto A2c = A2(omega);
    auto A = fd_op->GetSystemMatrix(1.0 + 0.0i, lambda, lambda * lambda, K.get(), C0,
                                    M.get(), A2c.get());
    auto P = fd_op->GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, lambda,
                                                             lambda * lambda, omega);
    const auto *mg = dynamic_cast<const ComplexMultigridOperator *>(P.get());
    REQUIRE(mg);
    A->Mult(x, y1);
    mg->GetFinestOperator().Mult(x, y2);
    CHECK_THAT(RelErr(y2, y1), Catch::Matchers::WithinAbs(0.0, 1.0e-12));
  }
}

TEST_CASE("SpaceOperator PML with unit stretch reproduces the bulk operators",
          "[spaceoperator][pml][Serial][Parallel]")
{
  using namespace std::complex_literals;

  // With σ = 0, κ = 1, the PML tensors equal the background material tensors, so all PML
  // terms (including the Floquet terms) must reproduce the standard bulk integrators.
  MPI_Comm comm = Mpi::World();
  mfem::Mesh serial_mesh = mfem::Mesh::MakeCartesian3D(3, 3, 3, mfem::Element::HEXAHEDRON);
  for (int i = 0; i < serial_mesh.GetNE(); i++)
  {
    mfem::Vector center;
    serial_mesh.GetElementCenter(i, center);
    serial_mesh.SetAttribute(i, (center(0) > 2.0 / 3.0 || center(2) > 2.0 / 3.0) ? 2 : 1);
  }
  serial_mesh.SetAttributes();
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(comm, serial_mesh));

  IntegrationSettingsGuard settings_guard;
  config::SolverData solver;
  solver.order = 2;
  solver.pa_order_threshold = 2;
  solver.linear.mg_max_levels = 2;
  solver.linear.mg_coarsening = MultigridCoarsening::LINEAR;
  solver.linear.pc_mat_shifted = 0;
  BilinearForm::pa_order_threshold = solver.pa_order_threshold;
  fem::DefaultIntegrationOrder::p_trial = solver.order;
  fem::DefaultIntegrationOrder::q_order_jac = solver.q_order_jac;
  fem::DefaultIntegrationOrder::q_order_extra_pk = solver.q_order_extra;
  fem::DefaultIntegrationOrder::q_order_extra_qk = solver.q_order_extra;

  const bool isotropic = GENERATE(false, true);
  const bool scaled = GENERATE(false, true);
  auto MakeSpaceOperator = [&](bool with_pml)
  {
    config::MaterialData vacuum, layer;
    vacuum.attributes = {1};
    layer.attributes = {2};
    if (isotropic)
    {
      layer.mu_r.s = {1.5, 1.5, 1.5};
      layer.epsilon_r.s = {3.0, 3.0, 3.0};
      layer.tandelta.s = {2.0e-2, 2.0e-2, 2.0e-2};
    }
    else
    {
      // Rotated material axes, so the background tensors are full.
      const double c = std::cos(0.3), s = std::sin(0.3);
      const std::array<std::array<double, 3>, 3> axes = {
          {{c, s, 0.0}, {-s * c, c * c, s}, {s * s, -c * s, c}}};
      layer.mu_r.s = {1.0, 2.0, 3.0};
      layer.mu_r.v = axes;
      layer.epsilon_r.s = {2.0, 4.0, 5.0};
      layer.epsilon_r.v = axes;
      layer.tandelta.s = {1.0e-2, 2.0e-2, 3.0e-2};
      layer.tandelta.v = axes;
    }
    if (with_pml)
    {
      layer.pml = config::PMLData();
      layer.pml->autodetect_geometry = false;
      layer.pml->direction_signs = {0, 1, 0, 0, 0, 1};
      layer.pml->thickness = {0.0, 1.0 / 3.0, 0.0, 0.0, 0.0, 1.0 / 3.0};
      layer.pml->sigma_max = {0.0, 0.0, 0.0};
      layer.pml->reference_frequency = 2.0;
    }
    config::DomainData domains;
    domains.attributes = {1, 2};
    domains.materials = {vacuum, layer};
    config::BoundaryData boundaries;
    boundaries.periodic.wave_vector = {0.3, -0.2, 0.5};
    boundaries.periodic.floquet_reference_freq = scaled ? 1.5 : 0.0;
    Units units(1.0, 1.0);
    return std::make_unique<SpaceOperator>(
        solver, domains, boundaries, scaled ? ProblemType::DRIVEN : ProblemType::EIGENMODE,
        units, mesh);
  };
  auto pml_op = MakeSpaceOperator(true), ref_op = MakeSpaceOperator(false);
  REQUIRE(pml_op->GetMaterialOp().HasPML());
  REQUIRE(pml_op->GetMaterialOp().HasWaveVector());
  REQUIRE(pml_op->GetMaterialOp().HasFloquetFrequencyScaling() == scaled);
  REQUIRE(!ref_op->GetMaterialOp().HasPML());

  auto CheckEqual = [&](const auto *A, const auto *B, auto &x, auto &y1, auto &y2)
  {
    REQUIRE(A);
    REQUIRE(B);
    for (int seed : {42, 314159})
    {
      if constexpr (std::is_same_v<std::decay_t<decltype(x)>, ComplexVector>)
      {
        x.Real().Randomize(seed + Mpi::Rank(comm));
        x.Imag().Randomize(2 * seed + 1 + Mpi::Rank(comm));
      }
      else
      {
        x.Randomize(seed + Mpi::Rank(comm));
      }
      A->Mult(x, y1);
      B->Mult(x, y2);
      const double norm = linalg::Norml2(comm, y2);
      REQUIRE(norm > 0.0);
      y1.Add(-1.0, y2);
      CHECK_THAT(linalg::Norml2(comm, y1) / norm, Catch::Matchers::WithinAbs(0.0, 1.0e-13));
    }
  };
  auto CheckComplex = [&](const ComplexOperator *A, const ComplexOperator *B)
  {
    ComplexVector x(A->Width()), y1(A->Height()), y2(A->Height());
    x.UseDevice(true);
    y1.UseDevice(true);
    y2.UseDevice(true);
    CheckEqual(A, B, x, y1, y2);
  };

  SECTION("System matrices")
  {
    constexpr auto diag = Operator::DIAG_ONE;
    CheckComplex(pml_op->GetStiffnessMatrix<ComplexOperator>(diag).get(),
                 ref_op->GetStiffnessMatrix<ComplexOperator>(diag).get());
    CheckComplex(pml_op->GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO).get(),
                 ref_op->GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO).get());
    if (scaled)
    {
      CheckComplex(pml_op->GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO).get(),
                   ref_op->GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO).get());
    }
    CHECK(
        !pml_op->GetExtraSystemMatrix(std::complex<double>(2.0, 0.1), Operator::DIAG_ZERO));
  }

  SECTION("Preconditioner")
  {
    // The real-valued approximation uses the entrywise magnitudes of the PML tensors, which
    // matches the bulk approximation for isotropic materials.
    for (bool pc_mat_real : {false, true})
    {
      if (pc_mat_real && !isotropic)
      {
        continue;
      }
      solver.linear.pc_mat_real = pc_mat_real;
      auto pml_pc_op = MakeSpaceOperator(true), ref_pc_op = MakeSpaceOperator(false);
      const std::complex<double> omega = {2.0, 0.1}, lambda = 1i * omega;
      auto P1 = pml_pc_op->GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, lambda,
                                                                    lambda * lambda, omega);
      auto P2 = ref_pc_op->GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, lambda,
                                                                    lambda * lambda, omega);
      const auto *mg1 = dynamic_cast<const ComplexMultigridOperator *>(P1.get());
      const auto *mg2 = dynamic_cast<const ComplexMultigridOperator *>(P2.get());
      REQUIRE(mg1);
      REQUIRE(mg2);
      REQUIRE(mg1->GetNumLevels() == 2);
      for (std::size_t l = 0; l < mg1->GetNumLevels(); l++)
      {
        CheckComplex(&mg1->GetOperatorAtLevel(l), &mg2->GetOperatorAtLevel(l));
        CheckComplex(&mg1->GetAuxiliaryOperatorAtLevel(l),
                     &mg2->GetAuxiliaryOperatorAtLevel(l));
      }
    }
  }
}

TEST_CASE("SpaceOperator PML operators match an MFEM coefficient reference",
          "[spaceoperator][pml][Serial][Parallel]")
{
  using namespace std::complex_literals;

  // Unit cube with a PML layer on the +x and +z faces (including the edge region), on
  // tetrahedral and hexahedral meshes. The PML curl-curl and mass operators are compared
  // to MFEM integrators with the PML tensors evaluated at the quadrature points by an
  // independent implementation of the stretch.
  const auto element = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  MPI_Comm comm = Mpi::World();
  auto serial_mesh = mfem::Mesh::MakeCartesian3D(3, 3, 3, element);
  for (int i = 0; i < serial_mesh.GetNE(); i++)
  {
    mfem::Vector center;
    serial_mesh.GetElementCenter(i, center);
    serial_mesh.SetAttribute(i, (center(0) > 2.0 / 3.0 || center(2) > 2.0 / 3.0) ? 2 : 1);
  }
  serial_mesh.SetAttributes();
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(comm, serial_mesh));

  IntegrationSettingsGuard settings_guard;
  config::SolverData solver;
  solver.order = 2;
  solver.linear.mg_max_levels = 1;
  fem::DefaultIntegrationOrder::p_trial = solver.order;

  config::MaterialData vacuum, layer;
  vacuum.attributes = {1};
  layer.attributes = {2};
  layer.mu_r.s = {1.0, 2.0, 1.5};
  layer.epsilon_r.s = {2.0, 3.0, 4.0};
  layer.tandelta.s = {1.0e-2, 2.0e-2, 3.0e-2};
  layer.pml = config::PMLData();
  layer.pml->autodetect_geometry = false;
  layer.pml->direction_signs = {0, 1, 0, 0, 0, 1};
  layer.pml->thickness = {0.0, 1.0 / 3.0, 0.0, 0.0, 0.0, 1.0 / 3.0};
  layer.pml->sigma_max = {3.0, 0.0, 5.0};
  layer.pml->kappa_max = {1.5, 1.0, 2.0};
  layer.pml->alpha_max = {0.2, 0.0, 0.1};
  layer.pml->reference_frequency = 2.0;
  config::DomainData domains;
  domains.attributes = {1, 2};
  domains.materials = {vacuum, layer};
  config::BoundaryData boundaries;
  Units units(1.0, 1.0);
  SpaceOperator space_op(solver, domains, boundaries, ProblemType::DRIVEN, units, mesh);
  const auto &profiles = space_op.GetMaterialOp().GetPMLProfiles();
  REQUIRE(profiles.size() == 1);
  const auto &p = profiles[0];

  // Reference stretch factors and tensors at a physical point (diagonal background).
  auto Stretch = [&p](const mfem::Vector &x)
  {
    std::array<std::complex<double>, 3> s;
    const std::array<double, 3> xp = {x(0), x(1), x(2)};
    const auto r = pml::ComputeDepthFraction(p, xp);
    for (int a = 0; a < 3; a++)
    {
      const double shape = std::pow(r[a], p.order);
      const bool pos = (xp[a] > p.geometry.inner[2 * a + 1]);
      const double sigma = shape * p.sigma_max[2 * a + (pos ? 1 : 0)];
      s[a] = 1.0 + (p.kappa_max[a] - 1.0) * shape +
             sigma / (p.alpha_max[a] * shape + 1i * p.reference_frequency);
    }
    return s;
  };
  auto MakeCoefficient = [&](bool muinv, bool imag)
  {
    return mfem::MatrixFunctionCoefficient(
        3,
        [=, &p](const mfem::Vector &x, mfem::DenseMatrix &T)
        {
          // Bulk material values are evaluated by the PWConstCoefficient below.
          const auto s = Stretch(x);
          const auto det = s[0] * s[1] * s[2];
          T = 0.0;
          for (int i = 0; i < 3; i++)
          {
            const std::complex<double> b =
                muinv ? p.mu_inv[4 * i]
                      : p.epsilon_real[4 * i] + 1i * p.epsilon_imag[4 * i];
            const std::complex<double> t =
                muinv ? b * s[i] * s[i] / det : b * det / (s[i] * s[i]);
            T(i, i) = imag ? t.imag() : t.real();
          }
        });
  };

  const auto &pfes = space_op.GetNDSpace().Get();
  const int q_order = fem::DefaultIntegrationOrder::Get(
      mesh[0]->Get(), mesh[0]->Get().GetElementGeometry(0));
  auto AssembleReference = [&](bool muinv, bool imag)
  {
    mfem::ParBilinearForm a(const_cast<mfem::ParFiniteElementSpace *>(&pfes));
    mfem::Array<int> pml_marker({0, 1}), bulk_marker({1, 0});
    auto *pml_coef = new mfem::MatrixFunctionCoefficient(MakeCoefficient(muinv, imag));
    auto *bulk_coef = new mfem::ConstantCoefficient(muinv ? 1.0 : 1.0);
    const auto &ir = mfem::IntRules.Get(mesh[0]->Get().GetElementGeometry(0), q_order);
    mfem::BilinearFormIntegrator *integ_pml, *integ_bulk = nullptr;
    if (muinv)
    {
      integ_pml = new mfem::CurlCurlIntegrator(*pml_coef);
      if (!imag)
      {
        integ_bulk = new mfem::CurlCurlIntegrator(*bulk_coef);
      }
    }
    else
    {
      integ_pml = new mfem::VectorFEMassIntegrator(*pml_coef);
      if (!imag)
      {
        integ_bulk = new mfem::VectorFEMassIntegrator(*bulk_coef);
      }
    }
    integ_pml->SetIntRule(&ir);
    a.AddDomainIntegrator(integ_pml, pml_marker);
    if (integ_bulk)
    {
      integ_bulk->SetIntRule(&ir);
      a.AddDomainIntegrator(integ_bulk, bulk_marker);
    }
    a.Assemble();
    a.Finalize();
    std::unique_ptr<mfem::HypreParMatrix> A(a.ParallelAssemble());
    delete pml_coef;
    delete bulk_coef;
    return A;
  };

  for (bool muinv : {true, false})
  {
    auto Ar = AssembleReference(muinv, false), Ai = AssembleReference(muinv, true);
    auto A = muinv ? space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO)
                   : space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    REQUIRE(A);
    ComplexVector x(A->Width()), y(A->Height()), y_ref(A->Height());
    x.Real().Randomize(42 + Mpi::Rank(comm));
    x.Imag().Randomize(4242 + Mpi::Rank(comm));
    A->Mult(x, y);
    Vector t(A->Height());
    Ar->Mult(x.Real(), y_ref.Real());
    Ai->Mult(x.Imag(), t);
    y_ref.Real().Add(-1.0, t);
    Ar->Mult(x.Imag(), y_ref.Imag());
    Ai->Mult(x.Real(), t);
    y_ref.Imag().Add(1.0, t);
    const double norm = linalg::Norml2(comm, y_ref);
    y.Add(-1.0, y_ref);
    INFO("muinv = " << muinv << ", element = " << static_cast<int>(element));
    CHECK_THAT(linalg::Norml2(comm, y) / norm, Catch::Matchers::WithinAbs(0.0, 1.0e-12));
  }
}

}  // namespace palace
