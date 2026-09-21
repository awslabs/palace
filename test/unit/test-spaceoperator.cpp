// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <complex>
#include <limits>
#include <typeinfo>
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

// The real and imaginary libCEED leaves of a partially assembled complex operator may be
// packed into one fused two-component application, which ceed::CreateComplexOperator
// returns as a private subclass of ComplexWrapperOperator: its dynamic type identifies it.
// The packing itself, and the backend conditions it needs, are covered by the ComplexPacked
// tests in test-libceed.cpp.
bool IsFusedComplexApplication(const ComplexOperator &op)
{
  return typeid(op) != typeid(ComplexWrapperOperator);
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
  CHECK(std::all_of(zero_data.coarse_values.begin(), zero_data.coarse_values.end(),
                    [](double value) { return value == 0.0; }));
  CHECK(std::any_of(tiny_data.coarse_values.begin(), tiny_data.coarse_values.end(),
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

// The driven sweep passes the preconditioner to the Krylov solver as the system operator
// when its finest level applies exactly A = K + iω C - ω² (Mr + i Mi) (see
// SpaceOperator::PreconditionerMatchesDrivenSystemMatrix, SpaceOperator::ApplySameMatrix
// and DrivenSolver::SweepUniform).
TEST_CASE("Driven preconditioner finest level is the system matrix",
          "[spaceoperator][Serial][Parallel]")
{
  using namespace std::complex_literals;

  const auto comm = Mpi::World();
  auto serial_mesh = mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::TETRAHEDRON);
  // A material loss tangent and an electrical conductivity cannot be combined on the same
  // attribute, so split the domain: the first material gives the mass term an imaginary
  // part, the second gives a damping term, and both are needed for every part of A to be
  // nonzero.
  for (int i = 0; i < serial_mesh.GetNE(); i++)
  {
    serial_mesh.SetAttribute(i, 1 + (i % 2));
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
  solver.linear.pc_mat_real = false;
  solver.linear.pc_mat_shifted = 0;
  BilinearForm::pa_order_threshold = solver.pa_order_threshold;
  fem::DefaultIntegrationOrder::p_trial = solver.order;
  fem::DefaultIntegrationOrder::q_order_jac = solver.q_order_jac;
  fem::DefaultIntegrationOrder::q_order_extra_pk = solver.q_order_extra;
  fem::DefaultIntegrationOrder::q_order_extra_qk = solver.q_order_extra;

  config::MaterialData lossy;
  lossy.attributes = {1};
  lossy.mu_r.s = {1.0, 2.0, 3.0};
  lossy.epsilon_r.s = {2.0, 4.0, 5.0};
  lossy.tandelta.s = {0.1, 0.1, 0.1};
  config::MaterialData conductive;
  conductive.attributes = {2};
  conductive.epsilon_r.s = {3.0, 3.0, 3.0};
  conductive.sigma.s = {0.5, 0.5, 0.5};
  config::DomainData domains;
  domains.attributes = {1, 2};
  domains.materials = {lossy, conductive};

  // Essential dofs, so the diagonal policies of the two operators are exercised, and a
  // surface impedance boundary contributing to the boundary stiffness, damping and mass.
  config::ImpedanceData impedance;
  impedance.Rs = 2.0;
  impedance.Ls = 0.5;
  impedance.Cs = 0.25;
  impedance.attributes = {4};
  config::BoundaryData boundaries;
  boundaries.pec.attributes = {1, 2, 3};
  boundaries.impedance = {impedance};
  Units units(1.0, 1.0);

  double omega = 2.0;
  // Assemble the system matrix from the fixed K, C and M and the preconditioner with the
  // same driven coefficients, and return the relative difference of their application to a
  // random complex vector. The Krylov solver applies the preconditioner through the
  // multigrid wrapper, which forwards to its finest level.
  auto ApplyDifference = [&](SpaceOperator &space_op)
  {
    auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
    auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    REQUIRE(K);
    REQUIRE(C);
    REQUIRE(M);
    auto A2 = space_op.GetExtraSystemOperator(omega, Operator::DIAG_ZERO);
    CHECK(!A2);
    auto A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i,
                                      K.get(), C.get(), M.get(), A2.get());
    auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(
        1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i, omega);
    const auto *mg = dynamic_cast<const ComplexMultigridOperator *>(P.get());
    REQUIRE(mg);
    REQUIRE(mg->GetNumLevels() == 2);
    REQUIRE(&mg->GetFinestOperator() == &mg->GetOperatorAtLevel(1));

    ComplexVector x(A->Width()), expected(A->Height()), result(A->Height());
    x.UseDevice(true);
    expected.UseDevice(true);
    result.UseDevice(true);
    double max_rel_diff = 0.0;
    for (int seed : {42, 314159})
    {
      x.Real().Randomize(seed + Mpi::Rank(comm));
      x.Imag().Randomize(seed + 1 + Mpi::Rank(comm));
      A->Mult(x, expected);
      P->Mult(x, result);
      CHECK(linalg::Norml2(comm, expected) > 0.0);
      result.AXPY(-1.0, expected);
      max_rel_diff = std::max(max_rel_diff, linalg::Norml2(comm, result) /
                                                linalg::Norml2(comm, expected));
    }
    // The shipped gate must agree with the comparison done here.
    CHECK(space_op.ApplySameMatrix(*A, *P) == (max_rel_diff <= 1.0e-11));
    return max_rel_diff;
  };

  // With the sharing enabled the Krylov solver drives one object from both roles within an
  // iteration: the system apply through the multigrid wrapper and the residual/smoother
  // applies of the V-cycle on the same finest level. That level owns the leaves its fused
  // two-component application packs, and only its forward Mult/AddMult use the packed
  // composite, the transpose using the retained originals. The preconditioner stays the
  // sole owner: both roles are non-owning views of it. Check that the shared object is a
  // full stand-in for the K + iω C - ω² (Mr + i Mi) sum in either role, and that the
  // mutable staging buffers of the fused application carry no state between applications.
  auto CheckSharedApplication = [&](SpaceOperator &space_op)
  {
    auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
    auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i,
                                      K.get(), C.get(), M.get(),
                                      static_cast<const ComplexOperator *>(nullptr));
    auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(
        1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i, omega);
    const auto *mg = dynamic_cast<const ComplexMultigridOperator *>(P.get());
    const auto *system = dynamic_cast<const ComplexParOperator *>(A.get());
    REQUIRE(mg);
    REQUIRE(system);
    const auto *finest = dynamic_cast<const ComplexParOperator *>(&mg->GetFinestOperator());
    REQUIRE(finest);
    // The system matrix is always the plain wrapper around the two SumOperators of the
    // fixed K, C and M leaves, which is the application the sharing replaces by the finest
    // level's single (fused, where the backend packs it) one.
    CHECK(!IsFusedComplexApplication(system->LocalOperator()));

    ComplexVector x(A->Width()), y_system(A->Height()), y_level(A->Height()),
        y_repeat(A->Height());
    x.UseDevice(true);
    y_system.UseDevice(true);
    y_level.UseDevice(true);
    y_repeat.UseDevice(true);
    linalg::SetRandom(comm, x);
    auto Difference = [&comm](ComplexVector &y, const ComplexVector &y_ref)
    {
      y.AXPY(-1.0, y_ref);
      return linalg::Norml2(comm, y);
    };

    // The system role (through the multigrid wrapper) and the preconditioner's own use of
    // the level reach the same application, in either order and repeatedly.
    P->Mult(x, y_system);
    finest->Mult(x, y_level);
    P->Mult(x, y_repeat);
    CHECK(linalg::Norml2(comm, y_system) > 0.0);
    CHECK(Difference(y_level, y_system) == 0.0);
    CHECK(Difference(y_repeat, y_system) == 0.0);

    // The transpose goes through the leaves the fused application retains, and is the
    // transpose of the same matrix.
    ComplexVector yt(A->Width()), yt_ref(A->Width());
    yt.UseDevice(true);
    yt_ref.UseDevice(true);
    A->MultTranspose(x, yt_ref);
    P->MultTranspose(x, yt);
    CHECK(linalg::Norml2(comm, yt_ref) > 0.0);
    CHECK_THAT(Difference(yt, yt_ref) / linalg::Norml2(comm, yt_ref),
               Catch::Matchers::WithinAbs(0.0, 1.0e-12));
  };

  SECTION("Preconditioner is the system matrix")
  {
    SpaceOperator space_op(solver, domains, boundaries, ProblemType::DRIVEN, units, mesh);
    REQUIRE(space_op.PreconditionerMatchesDrivenSystemMatrix());
    REQUIRE(!space_op.HasFrequencyDependentBoundaryTerms());
    CHECK_THAT(ApplyDifference(space_op), Catch::Matchers::WithinAbs(0.0, 1.0e-13));
    CheckSharedApplication(space_op);
  }

  SECTION("Multiple attributes in one impedance boundary group")
  {
    // Two terms (1/Ls and -ω² Cs) accumulate into the same boundary coefficient entry,
    // which the impedance operator adds to once per attribute: the preconditioner counts
    // the second term once per attribute of the shared entry and is not the system matrix
    // (see MaterialPropertyCoefficient::AddMaterialProperty). The configuration conditions
    // cannot see this, so the numerical check is what disables the reuse.
    auto multi_attr_boundaries = boundaries;
    multi_attr_boundaries.impedance.front().attributes = {4, 5, 6};
    SpaceOperator space_op(solver, domains, multi_attr_boundaries, ProblemType::DRIVEN,
                           units, mesh);
    CHECK(space_op.PreconditionerMatchesDrivenSystemMatrix());
    CHECK(ApplyDifference(space_op) > 1.0e-6);
  }

  SECTION("Zero frequency makes the numerical check vacuous")
  {
    // Every term the multi-attribute mis-assembly above can over-count scales with omega
    // (damping) or omega^2 (mass): the first stamping into each coefficient part is always
    // correct and only the later, frequency-scaled ones are over-counted. At omega = 0 the
    // preconditioner and the system matrix therefore agree exactly even for this
    // mismatched configuration, which is why the driven sweeps keep the check armed until
    // the first nonzero frequency.
    auto multi_attr_boundaries = boundaries;
    multi_attr_boundaries.impedance.front().attributes = {4, 5, 6};
    SpaceOperator space_op(solver, domains, multi_attr_boundaries, ProblemType::DRIVEN,
                           units, mesh);
    omega = 0.0;
    CHECK_THAT(ApplyDifference(space_op), Catch::Matchers::WithinAbs(0.0, 1.0e-13));
  }

  SECTION("Shifted preconditioner is a different matrix")
  {
    auto shifted_solver = solver;
    shifted_solver.linear.pc_mat_shifted = 1;
    SpaceOperator space_op(shifted_solver, domains, boundaries, ProblemType::DRIVEN, units,
                           mesh);
    CHECK(!space_op.PreconditionerMatchesDrivenSystemMatrix());
    CHECK(ApplyDifference(space_op) > 1.0e-6);
  }

  SECTION("Real-valued preconditioner")
  {
    auto real_solver = solver;
    real_solver.linear.pc_mat_real = true;
    SpaceOperator space_op(real_solver, domains, boundaries, ProblemType::DRIVEN, units,
                           mesh);
    CHECK(!space_op.PreconditionerMatchesDrivenSystemMatrix());
  }

  SECTION("Single multigrid level")
  {
    auto single_solver = solver;
    single_solver.linear.mg_max_levels = 1;
    SpaceOperator space_op(single_solver, domains, boundaries, ProblemType::DRIVEN, units,
                           mesh);
    REQUIRE(space_op.GetNDSpaces().GetNumLevels() == 1);
    CHECK(!space_op.PreconditionerMatchesDrivenSystemMatrix());
  }

  SECTION("Floquet wave vector")
  {
    auto floquet_boundaries = boundaries;
    floquet_boundaries.periodic.wave_vector = {0.1, 0.0, 0.0};
    SpaceOperator space_op(solver, domains, floquet_boundaries, ProblemType::DRIVEN, units,
                           mesh);
    REQUIRE(space_op.GetMaterialOp().HasWaveVector());
    CHECK(!space_op.PreconditionerMatchesDrivenSystemMatrix());
  }

  SECTION("Second-order farfield boundary is frequency dependent")
  {
    auto farfield_boundaries = boundaries;
    farfield_boundaries.impedance = {};
    farfield_boundaries.farfield.order = 2;
    farfield_boundaries.farfield.attributes = {4, 5, 6};
    SpaceOperator space_op(solver, domains, farfield_boundaries, ProblemType::DRIVEN, units,
                           mesh);
    CHECK(space_op.PreconditionerMatchesDrivenSystemMatrix());
    CHECK(space_op.HasFrequencyDependentBoundaryTerms());
    CHECK(space_op.GetExtraSystemOperator(omega, Operator::DIAG_ZERO));
  }
}

}  // namespace palace
