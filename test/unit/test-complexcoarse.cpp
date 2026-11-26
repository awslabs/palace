// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <limits>
#include <memory>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "linalg/complexcoarse.hpp"
#include "linalg/mumps.hpp"
#include "linalg/mumpscomplex.hpp"
#include "linalg/rap.hpp"
#include "linalg/solver.hpp"
#include "linalg/strumpack.hpp"
#include "linalg/strumpackcomplex.hpp"
#include "linalg/superlu.hpp"
#include "linalg/superlucomplex.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"

#if defined(MFEM_USE_STRUMPACK) || defined(MFEM_USE_SUPERLU) || defined(MFEM_USE_MUMPS)

namespace palace
{

namespace
{

// The sparse direct solvers of this build which have an implementation in complex
// arithmetic.
std::vector<LinearSolver> ComplexArithmeticSolvers()
{
  std::vector<LinearSolver> types;
#if defined(MFEM_USE_STRUMPACK)
  types.push_back(LinearSolver::STRUMPACK);
#endif
#if defined(MFEM_USE_SUPERLU)
  types.push_back(LinearSolver::SUPERLU);
#endif
#if defined(MFEM_USE_MUMPS)
  types.push_back(LinearSolver::MUMPS);
#endif
  return types;
}

// The solver of the real-valued approximation and the complex arithmetic solver of the
// measured choice with the given sparse direct solver, as in ConfigurePreconditionerSolver.
std::unique_ptr<Solver<ComplexOperator>> MakeRealSolver(LinearSolver type, MPI_Comm comm)
{
  std::unique_ptr<mfem::Solver> pc;
  const auto reorder = SymbolicFactorization::DEFAULT;
  switch (type)
  {
#if defined(MFEM_USE_STRUMPACK)
    case LinearSolver::STRUMPACK:
      pc = std::make_unique<StrumpackSolver>(comm, reorder, SparseCompression::NONE, 0.0, 0,
                                             0, true, 0);
      break;
#endif
#if defined(MFEM_USE_SUPERLU)
    case LinearSolver::SUPERLU:
      pc = std::make_unique<SuperLUSolver>(comm, reorder, false, true, 0);
      break;
#endif
#if defined(MFEM_USE_MUMPS)
    case LinearSolver::MUMPS:
      pc = std::make_unique<MumpsSolver>(comm, MatrixSymmetry::SYMMETRIC, reorder, 0.0,
                                         true, 0);
      break;
#endif
    default:
      MFEM_ABORT("Unexpected solver type!");
  }
  return std::make_unique<MfemWrapperSolver<ComplexOperator>>(std::move(pc), false, false,
                                                              false, true);
}

std::unique_ptr<Solver<ComplexOperator>> MakeComplexSolver(LinearSolver type, MPI_Comm comm)
{
  const auto reorder = SymbolicFactorization::DEFAULT;
  switch (type)
  {
#if defined(MFEM_USE_STRUMPACK)
    case LinearSolver::STRUMPACK:
      return std::make_unique<StrumpackComplexSolver>(
          comm, reorder, SparseCompression::NONE, 0.0, 0, 0, false, true, 0);
#endif
#if defined(MFEM_USE_SUPERLU)
    case LinearSolver::SUPERLU:
      return std::make_unique<SuperLUComplexSolver>(comm, reorder, false, false, true, 0);
#endif
#if defined(MFEM_USE_MUMPS)
    case LinearSolver::MUMPS:
      return std::make_unique<MumpsComplexSolver>(comm, MatrixSymmetry::SYMMETRIC, reorder,
                                                  0.0, false, true, 0);
#endif
    default:
      MFEM_ABORT("Unexpected solver type!");
  }
  return {};
}

// The measured choice between the real-valued approximation and the complex-valued matrix
// with the given sparse direct solver.
std::unique_ptr<ComplexCoarseSolver> MakeMeasuredSolver(LinearSolver type, MPI_Comm comm,
                                                        const FiniteElementSpace &fespace)
{
  return std::make_unique<ComplexCoarseSolver>(
      comm, [type, comm]() { return MakeRealSolver(type, comm); },
      [type, comm]() { return MakeComplexSolver(type, comm); }, false, MeshNoise(fespace));
}

// A real-valued solve whose result is not a number, as after a failed factorization.
class NotANumberSolver : public Solver<ComplexOperator>
{
public:
  void SetOperator(const ComplexOperator &op) override
  {
    height = op.Height();
    width = op.Width();
  }

  void Mult(const ComplexVector &x, ComplexVector &y) const override
  {
    y = std::numeric_limits<double>::quiet_NaN();
  }
};

// Nédélec space of the given order on a tetrahedral mesh of the unit cube, with the
// operator Ar + i Ai, Ar = k K + mr M, Ai = mi M (K curl-curl, M mass matrix, no essential
// boundary conditions), assembled as for the coarse level of geometric multigrid.
struct TestProblem
{
  Mesh mesh;
  mfem::ND_FECollection fec;
  FiniteElementSpace fespace;
  MaterialOperator mat_op;

  int order;

  TestProblem(MPI_Comm comm, mfem::Mesh &serial_mesh, int order = 1)
    : mesh(comm, serial_mesh), fec(order, 3), fespace(mesh, &fec),
      mat_op({MaterialAttributes()}, config::PeriodicBoundaryData(), ProblemType::DRIVEN,
             mesh),
      order(order)
  {
  }

  static config::MaterialData MaterialAttributes()
  {
    config::MaterialData material;
    material.attributes = {1};
    return material;
  }

  std::unique_ptr<Operator> AssemblePart(double curlcurl, double mass) const
  {
    fem::DefaultIntegrationOrder::p_trial = order;
    MaterialPropertyCoefficient df(mat_op.MaxCeedAttribute()), f(mat_op.MaxCeedAttribute());
    df.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetInvPermeability(),
                      curlcurl);
    f.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetPermittivityReal(), mass);
    BilinearForm a(fespace);
    if (curlcurl != 0.0)
    {
      a.AddDomainIntegrator<CurlCurlIntegrator>(df);
    }
    a.AddDomainIntegrator<VectorFEMassIntegrator>(f);
    return a.Assemble(false);
  }

  std::unique_ptr<ComplexOperator> Assemble(double k, double mr, double mi) const
  {
    return std::make_unique<ComplexParOperator>(AssemblePart(k, mr), AssemblePart(0.0, mi),
                                                fespace);
  }

  // An operator with an imaginary part only, i d M, as the correction which SpaceOperator
  // assembles for a coarse-level operator.
  std::unique_ptr<ComplexOperator> AssembleCorrection(double d) const
  {
    return std::make_unique<ComplexParOperator>(nullptr, AssemblePart(0.0, d), fespace);
  }
};

mfem::Mesh MakeTestMesh(int n)
{
  auto mesh = mfem::Mesh::MakeCartesian3D(n, n, n, mfem::Element::TETRAHEDRON);
  while (mesh.GetNE() < 4 * Mpi::Size(Mpi::World()))
  {
    mesh.UniformRefinement();
  }
  return mesh;
}

// Relative residual of x for A x = b.
double Residual(MPI_Comm comm, const ComplexOperator &A, const ComplexVector &x,
                const ComplexVector &b)
{
  ComplexVector r(b.Size());
  r.UseDevice(true);
  A.Mult(x, r);
  r -= b;
  return linalg::Norml2(comm, r) / linalg::Norml2(comm, b);
}

}  // namespace

TEST_CASE("ComplexCoarseSolver chooses from the measured error of the real-valued solve",
          "[complexcoarse][strumpack][superlu][mumps][Serial][Parallel]")
{
  // For Ar = M and Ai = c M, the real-valued solve B = (Ar + Ai)⁻¹ has the error
  // I - B A = c (1 - i) / (1 + c) I, so ‖Eᵏ w‖ / ‖w‖ = (√2 c / (1 + c))ᵏ: below 1 (the
  // real-valued approximation is kept) for c = 0.1, above 1 (the complex-valued matrix is
  // factored) for c = 10. For c = 0, the imaginary part is assembled as zeros.
  MPI_Comm comm = Mpi::World();
  const auto type = GENERATE(from_range(ComplexArithmeticSolvers()));
  auto serial_mesh = MakeTestMesh(3);
  TestProblem problem(comm, serial_mesh);
  const auto A_fine = problem.Assemble(0.0, 1.0, 0.1);
  const auto A_poor = problem.Assemble(0.0, 1.0, 10.0);
  const auto A_zero = problem.Assemble(0.0, 1.0, 0.0);
  auto Expected = [](double c, int k)
  { return std::pow(std::sqrt(2.0) * c / (1.0 + c), k); };

  const int n = A_fine->Height();
  ComplexVector b(n), x(n);
  b.UseDevice(true);
  x.UseDevice(true);
  b.Real().Randomize(1);
  b.Imag().Randomize(2);

  // An operator whose imaginary part is zero: the real part is factored, which is exact,
  // without a measurement.
  auto solver = MakeMeasuredSolver(type, comm, problem.fespace);
  solver->SetOperator(*A_zero);
  CHECK(solver->GetRealSolver());
  CHECK(solver->GetRealSolveError().empty());
  solver->Mult(b, x);
  CHECK(Residual(comm, *A_zero, x, b) < 1.0e-10);

  // The real-valued approximation is kept: three solves are measured, and it is applied.
  solver->SetOperator(*A_fine);
  CHECK(solver->GetRealSolver());
  CHECK(!solver->GetComplexSolver());
  REQUIRE(solver->GetRealSolveError().size() == 3);
  for (int k = 0; k < 3; k++)
  {
    CHECK(solver->GetRealSolveError()[k] ==
          Catch::Approx(Expected(0.1, k + 1)).epsilon(1.0e-8));
  }
  solver->Mult(b, x);
  CHECK(Residual(comm, *A_fine, x, b) == Catch::Approx(Expected(0.1, 1)).epsilon(1.0e-8));

  // A later operator of the same solver whose real-valued approximation is poor: the
  // measurement stops at the first error above 1, the real-valued solver is released and
  // the complex-valued matrix is factored.
  solver->SetOperator(*A_poor);
  CHECK(!solver->GetRealSolver());
  CHECK(solver->GetComplexSolver());
  REQUIRE(solver->GetRealSolveError().size() == 1);
  CHECK(solver->GetRealSolveError()[0] == Catch::Approx(Expected(10.0, 1)).epsilon(1.0e-8));
  solver->Mult(b, x);
  CHECK(Residual(comm, *A_poor, x, b) < 1.0e-10);

  // The choice is kept for the later operators of the solver (as for the later frequencies
  // of a sweep), which are not measured, also where the real-valued approximation would be
  // fine, and after an operator whose imaginary part is zero.
  for (const auto *A : {A_fine.get(), A_zero.get(), A_fine.get()})
  {
    solver->SetOperator(*A);
    CHECK((solver->GetComplexSolver() != nullptr) == (A != A_zero.get()));
    CHECK((solver->GetRealSolver() != nullptr) == (A == A_zero.get()));
    CHECK(solver->GetRealSolveError().size() == 1);
    solver->Mult(b, x);
    CHECK(Residual(comm, *A, x, b) < 1.0e-10);
  }

  // A new solver (as on the next mesh of an adaptive refinement) measures again, and its
  // complex arithmetic solver is never built.
  auto new_solver = MakeMeasuredSolver(type, comm, problem.fespace);
  new_solver->SetOperator(*A_fine);
  CHECK(!new_solver->GetComplexSolver());
  CHECK(new_solver->GetRealSolveError().size() == 3);
  new_solver->Mult(b, x);
  CHECK(Residual(comm, *A_fine, x, b) == Catch::Approx(Expected(0.1, 1)).epsilon(1.0e-8));
}

TEST_CASE("ComplexCoarseSolver with an operator for the real-valued approximation",
          "[complexcoarse][strumpack][superlu][mumps][Serial][Parallel]")
{
  // The operator given to the solver is the one assembled for the real-valued
  // approximation, M + 0.1 i M, and the complex-valued operator is M + c i M, the former
  // plus a correction (c - 0.1) i M. The real-valued solve B = (1.1 M)⁻¹ then has the error
  // I - B A = (0.1 - c i) / 1.1 I for the complex-valued operator: its magnitude is 0.203
  // for c = 0.2 (the real-valued approximation is kept) and 9.09 for c = 10 (the
  // complex-valued operator is factored, not the operator which was given).
  MPI_Comm comm = Mpi::World();
  const auto type = GENERATE(from_range(ComplexArithmeticSolvers()));
  auto serial_mesh = MakeTestMesh(3);
  TestProblem problem(comm, serial_mesh);
  const auto A_approx = problem.Assemble(0.0, 1.0, 0.1);
  auto Expected = [](double c, int k) { return std::pow(std::hypot(0.1, c) / 1.1, k); };

  const int n = A_approx->Height();
  ComplexVector b(n), x(n);
  b.UseDevice(true);
  x.UseDevice(true);
  b.Real().Randomize(1);
  b.Imag().Randomize(2);

  // The flag which the solver shares with the assembler of its operators.
  bool complex_chosen = false;
  auto solver = MakeMeasuredSolver(type, comm, problem.fespace);
  {
    const auto D = problem.AssembleCorrection(0.2 - 0.1);
    const auto A = problem.Assemble(0.0, 1.0, 0.2);
    solver->SetChoiceContext(&complex_chosen, D.get());
    solver->SetOperator(*A_approx);
    CHECK(!complex_chosen);
    CHECK(solver->GetRealSolver());
    CHECK(!solver->GetComplexSolver());
    REQUIRE(solver->GetRealSolveError().size() == 3);
    for (int k = 0; k < 3; k++)
    {
      CHECK(solver->GetRealSolveError()[k] ==
            Catch::Approx(Expected(0.2, k + 1)).epsilon(1.0e-8));
    }
    solver->Mult(b, x);
    CHECK(Residual(comm, *A, x, b) == Catch::Approx(Expected(0.2, 1)).epsilon(1.0e-8));
  }
  {
    const auto D = problem.AssembleCorrection(10.0 - 0.1);
    const auto A = problem.Assemble(0.0, 1.0, 10.0);
    solver->SetChoiceContext(&complex_chosen, D.get());
    solver->SetOperator(*A_approx);
    CHECK(complex_chosen);
    CHECK(!solver->GetRealSolver());
    CHECK(solver->GetComplexSolver());
    REQUIRE(solver->GetRealSolveError().size() == 1);
    CHECK(solver->GetRealSolveError()[0] ==
          Catch::Approx(Expected(10.0, 1)).epsilon(1.0e-8));
    solver->Mult(b, x);
    CHECK(Residual(comm, *A, x, b) < 1.0e-10);

    // The later operators are assembled for the complex-valued solve, without a
    // correction, and are not measured, by this solver and by a new solver whose operators
    // come from the same assembler.
    auto new_solver = MakeMeasuredSolver(type, comm, problem.fespace);
    for (auto *pc : {solver.get(), new_solver.get()})
    {
      pc->SetChoiceContext(&complex_chosen, nullptr);
      pc->SetOperator(*A);
      CHECK(!pc->GetRealSolver());
      CHECK(pc->GetComplexSolver());
      CHECK(pc->GetRealSolveError().size() == ((pc == solver.get()) ? 1 : 0));
      pc->Mult(b, x);
      CHECK(Residual(comm, *A, x, b) < 1.0e-10);
    }
  }
}

TEST_CASE("ComplexCoarseSolver switches when the measured error is not a number",
          "[complexcoarse][strumpack][superlu][mumps][Serial][Parallel]")
{
  // A failed real-valued factorization (for example of a singular Ar + Ai) gives a
  // measured error which is not a number: the complex-valued matrix is factored.
  MPI_Comm comm = Mpi::World();
  const auto type = ComplexArithmeticSolvers().front();
  auto serial_mesh = MakeTestMesh(3);
  TestProblem problem(comm, serial_mesh);
  const auto A = problem.Assemble(0.0, 1.0, 0.1);
  ComplexCoarseSolver solver(
      comm, []() { return std::make_unique<NotANumberSolver>(); }, [type, comm]()
      { return MakeComplexSolver(type, comm); }, false, MeshNoise(problem.fespace));
  solver.SetOperator(*A);
  REQUIRE(solver.GetRealSolveError().size() == 1);
  CHECK(std::isnan(solver.GetRealSolveError()[0]));
  CHECK(!solver.GetRealSolver());
  CHECK(solver.GetComplexSolver());

  const int n = A->Height();
  ComplexVector b(n), x(n);
  b.UseDevice(true);
  x.UseDevice(true);
  b.Real().Randomize(1);
  b.Imag().Randomize(2);
  solver.Mult(b, x);
  CHECK(Residual(comm, *A, x, b) < 1.0e-10);
}

TEST_CASE("ComplexCoarseSolver measurement does not depend on the number of processes",
          "[complexcoarse][strumpack][superlu][mumps][Serial][Parallel]")
{
  // The measurement vector (MeshNoise), the measured errors and the choice on the
  // processes of the run match those of the same problem on one process, also on a
  // periodic mesh (where the elements on the two sides of a periodic boundary see different
  // coordinates for the same degree of freedom). The operators are
  // Ar = K + M, Ai = c M: the real-valued solve damps the components with a large curl and
  // amplifies the gradients by √2 c / (1 + c) per solve, so the measure decays for c = 0.2
  // and exceeds 1 at the third solve for c = 50.
  const auto type = ComplexArithmeticSolvers().front();
  const bool periodic = GENERATE(false, true);
  const int order = GENERATE(1, 2);
  auto serial_mesh = MakeTestMesh(6);
  if (periodic)
  {
    serial_mesh.EnsureNodes();
    std::vector<mfem::Vector> translations = {mfem::Vector({1.0, 0.0, 0.0})};
    serial_mesh = mfem::Mesh::MakePeriodic(
        serial_mesh, serial_mesh.CreatePeriodicVertexMapping(translations));
  }
  mfem::VectorFunctionCoefficient g(3,
                                    [](const mfem::Vector &x, mfem::Vector &f)
                                    {
                                      f(0) = std::sin(3.0 * x(0) + x(1));
                                      f(1) = std::cos(2.0 * x(1) - x(2));
                                      f(2) = x(0) * x(2);
                                    });
  mfem::Vector zero_vec(3);
  zero_vec = 0.0;
  mfem::VectorConstantCoefficient zero(zero_vec);

  struct Measured
  {
    HYPRE_BigInt size;
    double sum, norm, dist;
    std::vector<double> err_fine, err_poor;
    bool complex_fine, complex_poor;
  };
  auto Measure = [&](MPI_Comm comm)
  {
    TestProblem problem(comm, serial_mesh, order);
    const auto w = MeshNoise(problem.fespace);
    Measured m;
    m.size = problem.fespace.GlobalTrueVSize();
    REQUIRE(w.Size() == problem.fespace.GetTrueVSize());
    int signs = 1;
    for (int i = 0; i < w.Size(); i++)
    {
      signs = signs && (w(i) == 1.0 || w(i) == -1.0);
    }
    Mpi::GlobalMin(1, &signs, comm);
    CHECK(signs == 1);
    m.sum = linalg::Sum(comm, w);
    mfem::ParGridFunction gf(&problem.fespace.Get());
    gf.SetFromTrueDofs(w);
    m.norm = gf.ComputeL2Error(zero);
    m.dist = gf.ComputeL2Error(g);
    for (const double c : {0.2, 50.0})
    {
      auto solver = MakeMeasuredSolver(type, comm, problem.fespace);
      solver->SetOperator(*problem.Assemble(1.0, 1.0, c));
      (c < 1.0 ? m.err_fine : m.err_poor) = solver->GetRealSolveError();
      (c < 1.0 ? m.complex_fine : m.complex_poor) = (solver->GetComplexSolver() != nullptr);
    }
    return m;
  };
  const auto par = Measure(Mpi::World());
  const auto ser = Measure(MPI_COMM_SELF);
  INFO("periodic " << periodic << ", order " << order);
  INFO("measured errors: c = 0.2: "
       << fmt::format("{:.6e}", fmt::join(par.err_fine, ", "))
       << "; c = 50: " << fmt::format("{:.6e}", fmt::join(par.err_poor, ", ")));
  CHECK(par.size == ser.size);
  // At order 1 (the coarsest level of p-multigrid) the vector is the same function for
  // every partition. At higher orders on tetrahedra, the basis of the degrees of freedom
  // of a face depends on the orientation of the face in the parallel mesh, so the vector
  // differs there and the measured errors agree only approximately.
  const double tol = (order == 1) ? 1.0e-8 : 1.0e-2;
  if (order == 1)
  {
    CHECK(par.sum == ser.sum);
    CHECK(par.norm == Catch::Approx(ser.norm).epsilon(1.0e-12));
    CHECK(par.dist == Catch::Approx(ser.dist).epsilon(1.0e-12));
  }
  REQUIRE(par.err_fine.size() == ser.err_fine.size());
  for (std::size_t k = 0; k < par.err_fine.size(); k++)
  {
    CHECK(par.err_fine[k] == Catch::Approx(ser.err_fine[k]).epsilon(tol));
  }
  REQUIRE(par.err_poor.size() == ser.err_poor.size());
  for (std::size_t k = 0; k < par.err_poor.size(); k++)
  {
    CHECK(par.err_poor[k] == Catch::Approx(ser.err_poor[k]).epsilon(tol));
  }
  CHECK(!par.complex_fine);
  CHECK(!ser.complex_fine);
  CHECK(par.complex_poor);
  CHECK(ser.complex_poor);
}

}  // namespace palace

#endif
