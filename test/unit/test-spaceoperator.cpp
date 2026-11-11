// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
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
  // With single-precision smoothing, the auxiliary space operators above the coarsest
  // level store their element matrices in single precision (on CPUs, for simplices).
  const bool single = GENERATE(false, true);
  CAPTURE(element, single);
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
  solver.linear.mg_smooth_single = single;
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
      // The main levels above the coarsest also get a single-precision approximation for
      // the Chebyshev search directions, while their own action stays exact.
      const bool single_aux = single && auxiliary && l > 0;
      const auto *approx =
          dynamic_cast<const ParOperator &>(actual).GetSinglePrecisionOperator();
      if (!single || auxiliary || l == 0)
      {
        CHECK(!approx);
      }
      else if (element == mfem::Element::TETRAHEDRON)
      {
        CHECK(approx);
      }
      Vector x(space.GetTrueVSize()), expected(x.Size()), result(x.Size());
      for (int seed : {42, 314159})
      {
        x.Randomize(seed + Mpi::Rank(comm));
        reference.Mult(x, expected);
        actual.Mult(x, result);
        result.Add(-1.0, expected);
        const double error = linalg::Norml2(comm, result) / linalg::Norml2(comm, expected);
        CAPTURE(auxiliary, l, error);
        CHECK(error < (single_aux ? 2.0e-6 : 2.0e-13));
        if (single_aux && element == mfem::Element::TETRAHEDRON)
        {
          CHECK(error > 1.0e-11);
        }
        if (approx)
        {
          approx->Mult(x, result);
          result.Add(-1.0, expected);
          const double approx_error =
              linalg::Norml2(comm, result) / linalg::Norml2(comm, expected);
          CAPTURE(approx_error);
          CHECK(approx_error > 1.0e-11);
          CHECK(approx_error < 2.0e-6);
        }
      }
    }
  }
}

// The driven solvers apply the system operator through the finest level of the
// preconditioner when that level is the sparse part of A = K + iω C - ω² (Mr + i Mi) +
// A2(ω) (see SpaceOperator::GetSystemOperator). The wave-port modal correction, which is
// the matrix-free part of A2, is covered in test-waveportoperator.cpp.
TEST_CASE("Frequency sweeps combine the diagonals of the preconditioner terms",
          "[spaceoperator][Serial][Parallel]")
{
  using namespace std::complex_literals;

  const auto comm = Mpi::World();
  auto serial_mesh = mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::TETRAHEDRON);
  // A lossy dielectric and a conductor give the preconditioner every term (stiffness,
  // damping and mass, with a complex mass), as does the surface impedance boundary.
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
  config::ImpedanceData impedance;
  impedance.Rs = 2.0;
  impedance.Ls = 0.5;
  impedance.Cs = 0.25;
  impedance.attributes = {4};
  config::BoundaryData boundaries;
  boundaries.pec.attributes = {1, 2, 3};
  boundaries.impedance = {impedance};
  Units units(1.0, 1.0);

  auto Assemble = [](SpaceOperator &space_op, double omega)
  {
    return space_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, 1i * omega,
                                                             -omega * omega + 0.0i, omega);
  };
  // Diagonals of the levels above the coarsest of the Nédélec and auxiliary hierarchies.
  auto Diagonals = [](const ComplexOperator &P)
  {
    const auto *mg = dynamic_cast<const ComplexMultigridOperator *>(&P);
    REQUIRE(mg);
    REQUIRE(mg->GetNumLevels() == 2);
    std::vector<ComplexVector> d;
    for (bool aux : {false, true})
    {
      const auto &op = aux ? mg->GetAuxiliaryOperatorAtLevel(1) : mg->GetOperatorAtLevel(1);
      d.emplace_back(op.Height());
      d.back().UseDevice(true);
      op.AssembleDiagonal(d.back());
    }
    return d;
  };

  // A sweep combines the diagonals of the terms from its second frequency on, while the
  // first preconditioner of a new operator assembles them: both must agree to round-off,
  // including the eliminated essential dofs.
  SpaceOperator sweep_op(solver, domains, boundaries, ProblemType::DRIVEN, units, mesh);
  for (double omega : {2.0, 3.0, 5.0})
  {
    const auto d = Diagonals(*Assemble(sweep_op, omega));
    SpaceOperator ref_op(solver, domains, boundaries, ProblemType::DRIVEN, units, mesh);
    const auto d_ref = Diagonals(*Assemble(ref_op, omega));
    for (std::size_t i = 0; i < d.size(); i++)
    {
      const double norm = linalg::Norml2(comm, d_ref[i]);
      CHECK(norm > 0.0);
      CHECK(linalg::Norml2(comm, d_ref[i].Imag()) > 0.0);
      ComplexVector diff(d[i]);
      diff.UseDevice(true);
      diff.AXPY(-1.0, d_ref[i]);
      CHECK_THAT(linalg::Norml2(comm, diff) / norm,
                 Catch::Matchers::WithinAbs(0.0, 1.0e-14));
    }
  }
}

TEST_CASE("Frequency sweeps combine the single-precision smoother element matrices of the "
          "preconditioner terms",
          "[spaceoperator][Serial][Parallel]")
{
  using namespace std::complex_literals;

  const auto comm = Mpi::World();
  auto serial_mesh = mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::TETRAHEDRON);
  // A lossy dielectric, a conductor and vacuum: the volume terms of the imaginary part are
  // over a subset of the elements.
  for (int i = 0; i < serial_mesh.GetNE(); i++)
  {
    serial_mesh.SetAttribute(i, 1 + (i % 3));
  }
  serial_mesh.SetAttributes();
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(comm, serial_mesh));

  // On affine elements, the curl-curl leaves of the stiffness term use a reduced basis with
  // fewer points than the curl-curl + mass leaves of the preconditioner (see
  // ceed::InitReducedBasis).
  IntegrationSettingsGuard settings_guard;
  config::SolverData solver;
  solver.order = GENERATE(2, 4);
  solver.pa_order_threshold = 2;
  solver.linear.mg_max_levels = 2;
  solver.linear.mg_coarsening = MultigridCoarsening::LINEAR;
  solver.linear.mg_smooth_single = true;
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
  config::MaterialData vacuum;
  vacuum.attributes = {3};
  config::DomainData domains;
  domains.attributes = {1, 2, 3};
  domains.materials = {lossy, conductive, vacuum};
  config::ImpedanceData impedance;
  impedance.Rs = 2.0;
  impedance.Ls = 0.5;
  impedance.Cs = 0.25;
  impedance.attributes = {4};
  config::BoundaryData boundaries;
  boundaries.pec.attributes = {1, 2, 3};
  boundaries.impedance = {impedance};
  Units units(1.0, 1.0);

  auto Assemble = [](SpaceOperator &space_op, double omega)
  {
    return space_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, 1i * omega,
                                                             -omega * omega + 0.0i, omega);
  };
  // Applications of the single-precision approximation of the Nédélec operator and of the
  // auxiliary operator (with single-precision element matrices) of the fine level, and of
  // the exact Nédélec operator, to the same random real vector: the real and imaginary
  // parts of the results are those of the actions of the real and imaginary parts of the
  // operators.
  auto Apply = [&](const ComplexOperator &P)
  {
    const auto *mg = dynamic_cast<const ComplexMultigridOperator *>(&P);
    REQUIRE(mg);
    REQUIRE(mg->GetNumLevels() == 2);
    const auto *nd = dynamic_cast<const ComplexParOperator *>(&mg->GetOperatorAtLevel(1));
    REQUIRE(nd);
    REQUIRE(nd->GetSinglePrecisionOperator());
    const auto &aux = mg->GetAuxiliaryOperatorAtLevel(1);
    std::vector<ComplexVector> y;
    for (const ComplexOperator *op :
         {nd->GetSinglePrecisionOperator(), &aux, static_cast<const ComplexOperator *>(nd)})
    {
      ComplexVector x(op->Width());
      x.UseDevice(true);
      x.Real().Randomize(42 + Mpi::Rank(comm));
      x.Imag() = 0.0;
      y.emplace_back(op->Height());
      y.back().UseDevice(true);
      op->Mult(x, y.back());
    }
    return y;
  };
  // Relative differences of the real and imaginary parts.
  auto RelativeDifferences = [&](const ComplexVector &y, const ComplexVector &y_ref)
  {
    std::array<double, 2> diff;
    for (int part = 0; part < 2; part++)
    {
      const auto &a = part ? y.Imag() : y.Real();
      const auto &a_ref = part ? y_ref.Imag() : y_ref.Real();
      const double norm = linalg::Norml2(comm, a_ref);
      REQUIRE(norm > 0.0);
      Vector d(a);
      d.UseDevice(true);
      d.Add(-1.0, a_ref);
      diff[part] = linalg::Norml2(comm, d) / norm;
    }
    return diff;
  };

  // A sweep forms the single-precision element matrices of both parts from those of the
  // terms from its second frequency on, the first preconditioner of a new operator
  // assembles them: both approximate the exact operator to single precision, and differ
  // from each other by the rounding of the terms.
  SpaceOperator sweep_op(solver, domains, boundaries, ProblemType::DRIVEN, units, mesh);
  bool first = true;
  for (double omega : {2.0, 3.0, 5.0})
  {
    const auto y = Apply(*Assemble(sweep_op, omega));
    SpaceOperator ref_op(solver, domains, boundaries, ProblemType::DRIVEN, units, mesh);
    const auto y_ref = Apply(*Assemble(ref_op, omega));
    for (int k = 0; k < 2; k++)
    {
      for (double diff : RelativeDifferences(y[k], y_ref[k]))
      {
        CHECK(diff < 1.0e-6);
        if (first)
        {
          CHECK(diff == 0.0);
        }
        else
        {
          CHECK(diff > 0.0);
        }
      }
    }
    for (double diff : RelativeDifferences(y[0], y_ref[2]))
    {
      CHECK(diff < 1.0e-6);
    }
    for (double diff : RelativeDifferences(y[2], y_ref[2]))
    {
      CHECK(diff == 0.0);
    }
    first = false;
  }
}

TEST_CASE("Driven system operator from the preconditioner",
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
  auto GetPreconditioner = [&](SpaceOperator &space_op)
  {
    return space_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, 1i * omega,
                                                             -omega * omega + 0.0i, omega);
  };

  // The system operator assembled from the fixed K, C and M, as the driven solvers do when
  // they cannot apply it through the preconditioner. A refers to the other operators.
  struct AssembledSystem
  {
    std::unique_ptr<ComplexOperator> K, C, M, A2, A;
  };
  auto AssembleSystem = [&](SpaceOperator &space_op, bool include_extra)
  {
    AssembledSystem sys;
    sys.K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
    sys.C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    sys.M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    REQUIRE(sys.K);
    REQUIRE(sys.C);
    REQUIRE(sys.M);
    if (include_extra)
    {
      sys.A2 = space_op.GetExtraSystemOperator(omega, Operator::DIAG_ZERO);
    }
    sys.A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i,
                                     sys.K.get(), sys.C.get(), sys.M.get(), sys.A2.get());
    return sys;
  };

  // Relative difference of the application of two operators, or of their transposes, to
  // random complex vectors.
  auto RelativeDifference = [&](const ComplexOperator &expected_op,
                                const ComplexOperator &actual_op, bool transpose = false)
  {
    ComplexVector x(expected_op.Width()), expected(expected_op.Height()),
        result(expected_op.Height());
    x.UseDevice(true);
    expected.UseDevice(true);
    result.UseDevice(true);
    double max_rel_diff = 0.0;
    for (int seed : {42, 314159})
    {
      x.Real().Randomize(seed + Mpi::Rank(comm));
      x.Imag().Randomize(seed + 1 + Mpi::Rank(comm));
      if (transpose)
      {
        expected_op.MultTranspose(x, expected);
        actual_op.MultTranspose(x, result);
      }
      else
      {
        expected_op.Mult(x, expected);
        actual_op.Mult(x, result);
      }
      const double norm = linalg::Norml2(comm, expected);
      CHECK(norm > 0.0);
      result.AXPY(-1.0, expected);
      max_rel_diff = std::max(max_rel_diff, linalg::Norml2(comm, result) / norm);
    }
    return max_rel_diff;
  };

  // The system operator from the preconditioner must match the assembled one. When there
  // are frequency-dependent boundary terms, also check that they change the system
  // operator, so that the match covers them.
  auto CheckSystemOperator = [&](SpaceOperator &space_op, bool has_extra)
  {
    auto sys = AssembleSystem(space_op, true);
    CHECK(static_cast<bool>(sys.A2) == has_extra);
    auto P = GetPreconditioner(space_op);
    auto A = space_op.GetSystemOperator(omega, *P);
    REQUIRE(A);
    CHECK_THAT(RelativeDifference(*sys.A, *A), Catch::Matchers::WithinAbs(0.0, 1.0e-12));
    CHECK_THAT(RelativeDifference(*sys.A, *A, true),
               Catch::Matchers::WithinAbs(0.0, 1.0e-12));
    if (has_extra)
    {
      auto sys_no_extra = AssembleSystem(space_op, false);
      CHECK(RelativeDifference(*sys.A, *sys_no_extra.A) > 1.0e-6);
    }
  };

  SECTION("Surface impedance boundary")
  {
    SpaceOperator space_op(solver, domains, boundaries, ProblemType::DRIVEN, units, mesh);
    CheckSystemOperator(space_op, false);
  }

  SECTION("Multiple attributes in one impedance boundary group")
  {
    // Repeated terms on one boundary group must update only the listed attributes even when
    // equal material properties initially share coefficient storage.
    auto multi_attr_boundaries = boundaries;
    multi_attr_boundaries.impedance.front().attributes = {4, 5, 6};
    SpaceOperator space_op(solver, domains, multi_attr_boundaries, ProblemType::DRIVEN,
                           units, mesh);
    CheckSystemOperator(space_op, false);
  }

  SECTION("Second-order farfield boundary")
  {
    auto farfield_boundaries = boundaries;
    farfield_boundaries.impedance = {};
    farfield_boundaries.farfield.order = 2;
    farfield_boundaries.farfield.attributes = {4, 5, 6};
    SpaceOperator space_op(solver, domains, farfield_boundaries, ProblemType::DRIVEN, units,
                           mesh);
    CheckSystemOperator(space_op, true);
  }

  SECTION("Surface conductivity boundary")
  {
    auto conductivity_boundaries = boundaries;
    auto &conductivity = conductivity_boundaries.conductivity.emplace_back();
    conductivity.sigma = 0.5;
    conductivity.attributes = {5};
    SpaceOperator space_op(solver, domains, conductivity_boundaries, ProblemType::DRIVEN,
                           units, mesh);
    CheckSystemOperator(space_op, true);
  }

  SECTION("Rational impedance boundary")
  {
    // Series RL impedance Z(s) = 0.5 s + 2.
    auto rational_boundaries = boundaries;
    auto &rational = rational_boundaries.rational_impedance.emplace_back();
    rational.num = {0.5, 2.0};
    rational.den = {1.0};
    rational.attributes = {6};
    SpaceOperator space_op(solver, domains, rational_boundaries, ProblemType::DRIVEN, units,
                           mesh);
    CheckSystemOperator(space_op, true);
  }

  SECTION("Shifted preconditioner is a different matrix")
  {
    auto shifted_solver = solver;
    shifted_solver.linear.pc_mat_shifted = 1;
    SpaceOperator space_op(shifted_solver, domains, boundaries, ProblemType::DRIVEN, units,
                           mesh);
    auto P = GetPreconditioner(space_op);
    CHECK(!space_op.GetSystemOperator(omega, *P));
    auto sys = AssembleSystem(space_op, true);
    CHECK(RelativeDifference(*sys.A, *P) > 1.0e-6);
  }

  SECTION("Real-valued preconditioner")
  {
    auto real_solver = solver;
    real_solver.linear.pc_mat_real = true;
    SpaceOperator space_op(real_solver, domains, boundaries, ProblemType::DRIVEN, units,
                           mesh);
    auto P = GetPreconditioner(space_op);
    CHECK(!space_op.GetSystemOperator(omega, *P));
  }

  SECTION("Single multigrid level")
  {
    // Its only preconditioner level is the assembled coarse operator; keep the system
    // application matrix-free.
    auto single_solver = solver;
    single_solver.linear.mg_max_levels = 1;
    SpaceOperator space_op(single_solver, domains, boundaries, ProblemType::DRIVEN, units,
                           mesh);
    REQUIRE(space_op.GetNDSpaces().GetNumLevels() == 1);
    auto P = GetPreconditioner(space_op);
    CHECK(!space_op.GetSystemOperator(omega, *P));
  }

  SECTION("Floquet wave vector")
  {
    // Periodic wave-vector terms use a separate assembly path which is outside the reuse
    // contract.
    auto floquet_boundaries = boundaries;
    floquet_boundaries.periodic.wave_vector = {0.1, 0.0, 0.0};
    SpaceOperator space_op(solver, domains, floquet_boundaries, ProblemType::DRIVEN, units,
                           mesh);
    REQUIRE(space_op.GetMaterialOp().HasWaveVector());
    auto P = GetPreconditioner(space_op);
    CHECK(!space_op.GetSystemOperator(omega, *P));
  }
}

}  // namespace palace
