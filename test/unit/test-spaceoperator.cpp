// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <complex>
#include <limits>
#include <map>
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

// The local rows of an assembled parallel matrix as maps from the global column index to
// the value, which is independent of the order of the entries within a row and of the local
// numbering of the off-diagonal columns.
using AssembledRows = std::vector<std::map<HYPRE_BigInt, double>>;

AssembledRows GetAssembledRows(const mfem::HypreParMatrix &A)
{
  hypre_CSRMatrix *csr = hypre_MergeDiagAndOffd((hypre_ParCSRMatrix *)A);
  hypre_CSRMatrixMigrate(csr, HYPRE_MEMORY_HOST);
  const HYPRE_Int *I = hypre_CSRMatrixI(csr);
  const HYPRE_BigInt *J = hypre_CSRMatrixBigJ(csr);
  const double *data = hypre_CSRMatrixData(csr);
  AssembledRows rows(hypre_CSRMatrixNumRows(csr));
  for (std::size_t i = 0; i < rows.size(); i++)
  {
    for (HYPRE_Int j = I[i]; j < I[i + 1]; j++)
    {
      rows[i][J[j]] = data[j];
    }
  }
  hypre_CSRMatrixDestroy(csr);
  return rows;
}

// The sparsity pattern of the local rows, as the global column indices of each of them.
std::vector<std::vector<HYPRE_BigInt>> GetPattern(const AssembledRows &rows)
{
  std::vector<std::vector<HYPRE_BigInt>> pattern(rows.size());
  for (std::size_t i = 0; i < rows.size(); i++)
  {
    for (const auto &[column, value] : rows[i])
    {
      pattern[i].push_back(column);
    }
  }
  return pattern;
}

// The assembled parallel matrix of one part of a complex parallel operator.
const mfem::HypreParMatrix &GetAssembledPart(const ComplexParOperator &op, bool imag)
{
  const auto *part = dynamic_cast<const ParOperator *>(imag ? op.Imag() : op.Real());
  REQUIRE(part);
  return part->ParallelAssemble();
}

struct ImaginaryHierarchyData
{
  AssembledRows coarse_rows;
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

  // The coarsest level is either assembled here or combined from the cached
  // frequency-independent term matrices, so compare it after the parallel assembly.
  ImaginaryHierarchyData data;
  constexpr bool imag = true;
  data.coarse_rows = GetAssembledRows(GetAssembledPart(*coarse_op, imag));

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

  auto AllValues = [](const AssembledRows &rows)
  {
    std::vector<double> values;
    for (const auto &row : rows)
    {
      for (const auto &[column, value] : row)
      {
        values.push_back(value);
      }
    }
    return values;
  };
  const auto zero_values = AllValues(zero_data.coarse_rows);
  const auto tiny_values = AllValues(tiny_data.coarse_rows);
  CHECK(GetPattern(zero_data.coarse_rows) == GetPattern(tiny_data.coarse_rows));
  CHECK(std::all_of(zero_values.begin(), zero_values.end(),
                    [](double value) { return value == 0.0; }));
  CHECK(std::any_of(tiny_values.begin(), tiny_values.end(),
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

// The Chebyshev smoothers ask every preconditioner level above the coarsest for its
// diagonal at each frequency, which the level operators assemble with libCEED. The diagonal
// is linear in the material property coefficients, so it is instead combined from the
// cached diagonals of the four frequency-independent terms (stiffness, damping, real and
// imaginary mass, with their boundary counterparts) using the same scalings as the
// coefficients themselves (SpaceOperator::GetPreconditionerScalars). A configuration whose
// coefficients are not those terms scaled by real numbers keeps the operator assembly.
TEST_CASE("Complex preconditioner diagonal is combined from cached terms",
          "[spaceoperator][Serial][Parallel]")
{
  using namespace std::complex_literals;

  const auto comm = Mpi::World();
  auto serial_mesh = mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::TETRAHEDRON);
  // A material loss tangent and an electrical conductivity cannot be combined on the same
  // attribute, so split the domain: both the damping and the imaginary mass term are then
  // nonzero, on different attributes.
  for (int i = 0; i < serial_mesh.GetNE(); i++)
  {
    serial_mesh.SetAttribute(i, 1 + (i % 2));
  }
  serial_mesh.SetAttributes();
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(comm, serial_mesh));

  IntegrationSettingsGuard settings_guard;
  config::SolverData solver;
  solver.order = 3;
  solver.pa_order_threshold = 2;
  solver.linear.mg_max_levels = 3;
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

  // Essential dofs, so the eliminated diagonal entries are exercised, and a surface
  // impedance boundary contributing to the boundary stiffness, damping and mass terms.
  config::ImpedanceData impedance;
  impedance.Rs = 2.0;
  impedance.Ls = 0.5;
  impedance.Cs = 0.25;
  impedance.attributes = {4};
  config::BoundaryData boundaries;
  boundaries.pec.attributes = {1, 2, 3};
  boundaries.impedance = {impedance};
  Units units(1.0, 1.0);

  constexpr double omega = 2.0;
  // Compare the diagonal which the smoothers get on every level above the coarsest of both
  // hierarchies with the one assembled from the level operator itself, and return the
  // largest relative infinity-norm difference of the two parts, which is at rounding level
  // (the largest measured here is 7.6e-16, a few ULP: the combination applies the term
  // scalings after the quadrature summation and the |P|ᵀ parallel assembly instead of
  // before). The type of a3 selects the
  // A2 stamping path exactly as for the solvers (double for driven, std::complex for the
  // eigenmode nonlinear solve).
  auto DiagonalDifference =
      [&](SpaceOperator &space_op, std::complex<double> lambda, auto a3, bool expect_cached)
  {
    auto pc = space_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, lambda,
                                                                lambda * lambda, a3);
    const auto *mg = dynamic_cast<const ComplexMultigridOperator *>(pc.get());
    REQUIRE(mg);
    REQUIRE(mg->GetNumLevels() == 3);
    double max_rel_diff = 0.0;
    for (bool aux : {false, true})
    {
      const auto &dbc_tdof_lists =
          aux ? space_op.GetH1DbcTDofLists() : space_op.GetNDDbcTDofLists();
      for (std::size_t l = 1; l < mg->GetNumLevels(); l++)
      {
        const auto *level_op = dynamic_cast<const ComplexParOperator *>(
            aux ? &mg->GetAuxiliaryOperatorAtLevel(l) : &mg->GetOperatorAtLevel(l));
        REQUIRE(level_op);
        CHECK(level_op->HasAssembledDiagonal() == expect_cached);

        ComplexVector diag(level_op->Height()), diag_ref(level_op->Height());
        diag.UseDevice(true);
        diag_ref.UseDevice(true);
        level_op->AssembleDiagonal(diag);
        level_op->AssembleDiagonalFromOperator(diag_ref);
        for (bool imag : {false, true})
        {
          const Vector &ref = imag ? diag_ref.Imag() : diag_ref.Real();
          Vector diff(imag ? diag.Imag() : diag.Real());
          diff -= ref;
          double norm[2] = {diff.Normlinf(), ref.Normlinf()};
          Mpi::GlobalMax(2, norm, comm);
          CHECK(norm[1] > 0.0);
          max_rel_diff = std::max(max_rel_diff, norm[0] / norm[1]);
        }

        // The essential dofs carry the eliminated values of the DIAG_ONE policy exactly.
        const auto &dbc_tdof_list = dbc_tdof_lists[l];
        int num_dbc = dbc_tdof_list.Size();
        Mpi::GlobalSum(1, &num_dbc, comm);
        REQUIRE(num_dbc > 0);
        const auto *dbc = dbc_tdof_list.HostRead();
        const auto *diag_real = diag.Real().HostRead();
        const auto *diag_imag = diag.Imag().HostRead();
        bool eliminated = true;
        for (int i = 0; i < dbc_tdof_list.Size(); i++)
        {
          eliminated = eliminated && diag_real[dbc[i]] == 1.0 && diag_imag[dbc[i]] == 0.0;
        }
        CHECK(eliminated);
      }
    }
    return max_rel_diff;
  };

  SECTION("Cached term diagonals")
  {
    // The driven coefficients at two frequencies: the second one reuses the cached terms.
    SpaceOperator space_op(solver, domains, boundaries, ProblemType::DRIVEN, units, mesh);
    constexpr bool expect_cached = true;
    CHECK(DiagonalDifference(space_op, 1i * omega, omega, expect_cached) <= 1.0e-13);
    CHECK(DiagonalDifference(space_op, 1.5i * omega, 1.5 * omega, expect_cached) <=
          1.0e-13);
  }

  SECTION("Complex frequency")
  {
    // A genuinely complex λ (the eigenmode nonlinear solve) scales every term into both
    // parts: a2 = λ² has an imaginary part, so the real mass term contributes to the
    // imaginary part and the imaginary mass term to the real one.
    const std::complex<double> lambda = 0.1 + 1i * omega;
    SpaceOperator space_op(solver, domains, boundaries, ProblemType::EIGENMODE, units,
                           mesh);
    constexpr bool expect_cached = true;
    CHECK(DiagonalDifference(space_op, lambda, lambda / 1i, expect_cached) <= 1.0e-13);
  }

  SECTION("Shifted preconditioner mass coefficient")
  {
    // The real mass term is scaled by |Re(a2)| instead of Re(a2), on both sides of the
    // shared scalings.
    auto shifted_solver = solver;
    shifted_solver.linear.pc_mat_shifted = 1;
    SpaceOperator space_op(shifted_solver, domains, boundaries, ProblemType::DRIVEN, units,
                           mesh);
    constexpr bool expect_cached = true;
    CHECK(DiagonalDifference(space_op, 1i * omega, omega, expect_cached) <= 1.0e-13);
  }

  SECTION("Floquet wave vector keeps the assembled diagonal")
  {
    // The periodic terms are stamped with their own frequency-dependent scalings.
    auto floquet_boundaries = boundaries;
    floquet_boundaries.periodic.wave_vector = {0.1, 0.0, 0.0};
    SpaceOperator space_op(solver, domains, floquet_boundaries, ProblemType::DRIVEN, units,
                           mesh);
    REQUIRE(space_op.GetMaterialOp().HasWaveVector());
    constexpr bool expect_cached = false;
    CHECK(DiagonalDifference(space_op, 1i * omega, omega, expect_cached) == 0.0);
  }

  SECTION("Frequency-dependent boundary keeps the assembled diagonal")
  {
    // The second-order farfield coefficients depend on ω non-linearly.
    auto farfield_boundaries = boundaries;
    farfield_boundaries.farfield.order = 2;
    farfield_boundaries.farfield.attributes = {5, 6};
    SpaceOperator space_op(solver, domains, farfield_boundaries, ProblemType::DRIVEN, units,
                           mesh);
    REQUIRE(space_op.HasFrequencyDependentBoundaryTerms());
    constexpr bool expect_cached = false;
    CHECK(DiagonalDifference(space_op, 1i * omega, omega, expect_cached) == 0.0);
  }
}

// The coarsest level of the complex preconditioner hierarchy is a sparse matrix which is
// linear in the frequency scalings of the four preconditioner terms (stiffness, damping,
// real and imaginary mass, with their boundary counterparts), so it is combined from the
// cached parallel matrices of those terms instead of being assembled and RAP'd at every
// frequency (SpaceOperator::CombinePreconditionerTermMatrices). The combination must have
// the same sparsity pattern as the per-frequency assembly, including the exact zeros of the
// terms which the current frequency does not scale in, the same values up to rounding, and
// the same essential true dof elimination.
TEST_CASE("Complex preconditioner coarsest level is combined from cached terms",
          "[spaceoperator][Serial][Parallel]")
{
  using namespace std::complex_literals;

  const auto comm = Mpi::World();
  auto serial_mesh = mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::TETRAHEDRON);
  // A material loss tangent and an electrical conductivity cannot be combined on the same
  // attribute, so split the domain: both the damping and the imaginary mass term are then
  // nonzero, on different attributes.
  for (int i = 0; i < serial_mesh.GetNE(); i++)
  {
    serial_mesh.SetAttribute(i, 1 + (i % 2));
  }
  serial_mesh.SetAttributes();
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(comm, serial_mesh));

  IntegrationSettingsGuard settings_guard;
  config::SolverData solver;
  solver.order = 3;
  solver.pa_order_threshold = 2;
  solver.linear.mg_max_levels = 3;
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

  // Essential dofs, so the eliminated rows are exercised, and a surface impedance boundary
  // contributing to the boundary stiffness, damping and mass terms.
  config::ImpedanceData impedance;
  impedance.Rs = 2.0;
  impedance.Ls = 0.5;
  impedance.Cs = 0.25;
  impedance.attributes = {4};
  config::BoundaryData boundaries;
  boundaries.pec.attributes = {1, 2, 3};
  boundaries.impedance = {impedance};
  Units units(1.0, 1.0);

  // The reference: a second-order farfield boundary on no attributes contributes to no
  // coefficient, but excludes the cached combination
  // (SpaceOperator::HasFrequencyDependentBoundaryTerms), so this operator assembles and
  // RAP's its coarsest level at every frequency, which is what the combination replaces.
  auto reference_boundaries = boundaries;
  reference_boundaries.farfield.order = 2;

  auto CoarseOperator = [](const ComplexOperator &pc)
  {
    const auto *mg = dynamic_cast<const ComplexMultigridOperator *>(&pc);
    REQUIRE(mg);
    REQUIRE(mg->GetNumLevels() == 3);
    const auto *coarse_op =
        dynamic_cast<const ComplexParOperator *>(&mg->GetOperatorAtLevel(0));
    REQUIRE(coarse_op);
    return coarse_op;
  };

  // Compare both parts of the combined coarsest level with the per-frequency assembly of
  // the reference, at the driven coefficients of one frequency.
  auto CheckCoarseLevel =
      [&](SpaceOperator &space_op, SpaceOperator &reference_op, double omega)
  {
    const std::complex<double> a1 = 1i * omega;
    auto pc =
        space_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, a1, a1 * a1, omega);
    auto pc_ref = reference_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, a1,
                                                                        a1 * a1, omega);
    const auto *coarse_op = CoarseOperator(*pc);
    const auto *coarse_ref_op = CoarseOperator(*pc_ref);
    CHECK(coarse_op->IsParallelAssembled());
    CHECK(!coarse_ref_op->IsParallelAssembled());

    const auto &dbc_tdof_list = space_op.GetNDDbcTDofLists()[0];
    int num_dbc = dbc_tdof_list.Size();
    Mpi::GlobalSum(1, &num_dbc, comm);
    REQUIRE(num_dbc > 0);
    double max_rel_diff = 0.0;
    for (bool imag : {false, true})
    {
      const auto &A = GetAssembledPart(*coarse_op, imag);
      const auto &A_ref = GetAssembledPart(*coarse_ref_op, imag);
      const auto rows = GetAssembledRows(A);
      const auto rows_ref = GetAssembledRows(A_ref);
      REQUIRE(rows.size() == rows_ref.size());
      CHECK(GetPattern(rows) == GetPattern(rows_ref));

      // The entries agree to rounding, relative to the largest entry of the row.
      for (std::size_t i = 0; i < rows.size(); i++)
      {
        double row_scale = 0.0, row_diff = 0.0;
        for (const auto &[column, value] : rows_ref[i])
        {
          row_scale = std::max(row_scale, std::abs(value));
          const auto it = rows[i].find(column);
          row_diff = std::max(row_diff,
                              std::abs((it != rows[i].end() ? it->second : 0.0) - value));
        }
        max_rel_diff =
            std::max(max_rel_diff, (row_scale > 0.0) ? row_diff / row_scale : row_diff);
      }

      // The essential dofs carry the eliminated values of the policy of the part exactly:
      // the real part uses DIAG_ONE and the imaginary part DIAG_ZERO.
      const auto *dbc = dbc_tdof_list.HostRead();
      bool eliminated = true;
      for (int i = 0; i < dbc_tdof_list.Size(); i++)
      {
        for (const auto &[column, value] : rows[dbc[i]])
        {
          const bool diagonal = (column == A.RowPart()[0] + dbc[i]);
          eliminated = eliminated && value == (diagonal && !imag ? 1.0 : 0.0);
        }
      }
      CHECK(eliminated);

      // The first entry in each row of the diagonal block is the diagonal one, as the
      // per-frequency assembly leaves it (ParOperator::ParallelAssemble).
      hypre_CSRMatrix *diag = hypre_ParCSRMatrixDiag((hypre_ParCSRMatrix *)A);
      hypre_CSRMatrixMigrate(diag, HYPRE_MEMORY_HOST);
      const HYPRE_Int *I = hypre_CSRMatrixI(diag);
      const HYPRE_Int *J = hypre_CSRMatrixJ(diag);
      bool diagonal_first = true;
      for (HYPRE_Int i = 0; i < hypre_CSRMatrixNumRows(diag); i++)
      {
        diagonal_first = diagonal_first && (I[i] == I[i + 1] || J[I[i]] == i);
      }
      CHECK(diagonal_first);
    }
    Mpi::GlobalMax(1, &max_rel_diff, comm);
    // Rounding level: the combination applies the term scalings after the assembly of the
    // term matrices instead of before (the largest measured here is 1.2e-15, a few ULP).
    CHECK(max_rel_diff <= 1.0e-13);

    // The combined level applies as the complex operator of its two assembled parts, which
    // is what the reference level applies matrix-free.
    ComplexVector x(coarse_op->Width()), y(coarse_op->Height()), y_ref(coarse_op->Height());
    x.UseDevice(true);
    y.UseDevice(true);
    y_ref.UseDevice(true);
    linalg::SetRandom(comm, x);
    coarse_op->Mult(x, y);
    coarse_ref_op->Mult(x, y_ref);
    y.AXPY(-1.0, y_ref);
    CHECK_THAT(linalg::Norml2(comm, y) / linalg::Norml2(comm, y_ref),
               Catch::Matchers::WithinAbs(0.0, 1.0e-12));
  };

  SECTION("Lossy materials")
  {
    // Every term has a coefficient, so every one of them has a cached matrix. The second
    // frequency reuses those matrices.
    SpaceOperator space_op(solver, domains, boundaries, ProblemType::DRIVEN, units, mesh);
    SpaceOperator reference_op(solver, domains, reference_boundaries, ProblemType::DRIVEN,
                               units, mesh);
    REQUIRE(!space_op.HasFrequencyDependentBoundaryTerms());
    REQUIRE(reference_op.HasFrequencyDependentBoundaryTerms());
    CheckCoarseLevel(space_op, reference_op, 2.0);
    CheckCoarseLevel(space_op, reference_op, 3.5);
  }

  SECTION("Coarse solver handoff")
  {
    // The coarse solver assembles both parts of the level and takes their matrices
    // (MfemWrapperSolver: ParallelAssemble, then StealParallelAssemble of each part and a
    // hypre sum of the two). Nothing is left of the level afterwards.
    SpaceOperator space_op(solver, domains, boundaries, ProblemType::DRIVEN, units, mesh);
    const std::complex<double> a1 = 2.0i;
    auto pc =
        space_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, a1, a1 * a1, 2.0);
    const auto *coarse_op = CoarseOperator(*pc);
    REQUIRE(coarse_op->IsParallelAssembled());
    const auto *real_part = dynamic_cast<const ParOperator *>(coarse_op->Real());
    const auto *imag_part = dynamic_cast<const ParOperator *>(coarse_op->Imag());
    REQUIRE(real_part);
    REQUIRE(imag_part);
    auto stolen_r = real_part->StealParallelAssemble();
    auto stolen_i = imag_part->StealParallelAssemble();
    REQUIRE(stolen_r);
    REQUIRE(stolen_i);
    CHECK(stolen_r->Height() == coarse_op->Height());
    std::unique_ptr<mfem::HypreParMatrix> sum(mfem::Add(1.0, *stolen_r, 1.0, *stolen_i));
    REQUIRE(sum);
    CHECK(sum->Height() == coarse_op->Height());
    CHECK_THROWS(real_part->ParallelAssemble());
  }

  SECTION("Lossless materials")
  {
    // Without a loss tangent or a conductivity the imaginary mass term has an exactly zero
    // coefficient and no cached matrix, while the per-frequency assembly of the level keeps
    // its configured, exactly zero block: the pattern must be the same either way.
    auto lossless = lossy;
    lossless.attributes = {1, 2};
    lossless.tandelta.s = {0.0, 0.0, 0.0};
    auto lossless_domains = domains;
    lossless_domains.materials = {lossless};
    SpaceOperator space_op(solver, lossless_domains, boundaries, ProblemType::DRIVEN, units,
                           mesh);
    SpaceOperator reference_op(solver, lossless_domains, reference_boundaries,
                               ProblemType::DRIVEN, units, mesh);
    REQUIRE(!space_op.GetMaterialOp().HasLossTangent());
    REQUIRE(!space_op.GetMaterialOp().HasConductivity());
    CheckCoarseLevel(space_op, reference_op, 2.0);
  }

  SECTION("Floquet wave vector keeps the assembled level")
  {
    // The periodic terms are stamped with their own frequency-dependent scalings, so this
    // level is assembled at every frequency as before.
    auto floquet_boundaries = boundaries;
    floquet_boundaries.periodic.wave_vector = {0.1, 0.0, 0.0};
    SpaceOperator space_op(solver, domains, floquet_boundaries, ProblemType::DRIVEN, units,
                           mesh);
    REQUIRE(space_op.GetMaterialOp().HasWaveVector());
    const std::complex<double> a1 = 2.0i;
    auto pc =
        space_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, a1, a1 * a1, 2.0);
    CHECK(!CoarseOperator(*pc)->IsParallelAssembled());
  }
}

}  // namespace palace
