// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <memory>
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

// Forget the measured choice of the complex-valued coarse solve before and after a test.
struct ComplexCoarseChoiceReset
{
  ComplexCoarseChoiceReset() { MeasuredComplexCoarseSolver::ResetChoice(); }
  ~ComplexCoarseChoiceReset() { MeasuredComplexCoarseSolver::ResetChoice(); }
};

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

// The measured choice between the real-valued approximation and the complex-valued matrix
// with the given sparse direct solver, as in ConfigurePreconditionerSolver.
std::unique_ptr<MeasuredComplexCoarseSolver>
MakeMeasuredSolver(LinearSolver type, MPI_Comm comm, const FiniteElementSpace &fespace)
{
  std::unique_ptr<mfem::Solver> real_pc;
  MeasuredComplexCoarseSolver::ComplexSolverFactory make_complex_pc;
  const auto reorder = SymbolicFactorization::DEFAULT;
  switch (type)
  {
#if defined(MFEM_USE_STRUMPACK)
    case LinearSolver::STRUMPACK:
      real_pc = std::make_unique<StrumpackSolver>(comm, reorder, SparseCompression::NONE,
                                                  0.0, 0, 0, true, 0);
      make_complex_pc = [comm, reorder]()
      { return std::make_unique<StrumpackComplexSolver>(comm, reorder, true, 0); };
      break;
#endif
#if defined(MFEM_USE_SUPERLU)
    case LinearSolver::SUPERLU:
      real_pc = std::make_unique<SuperLUSolver>(comm, reorder, false, true, 0);
      make_complex_pc = [comm, reorder]()
      { return std::make_unique<SuperLUComplexSolver>(comm, reorder, true, 0); };
      break;
#endif
#if defined(MFEM_USE_MUMPS)
    case LinearSolver::MUMPS:
      real_pc = std::make_unique<MumpsSolver>(comm, MatrixSymmetry::SYMMETRIC, reorder, 0.0,
                                              true, 0);
      make_complex_pc = [comm, reorder]()
      {
        return std::make_unique<MumpsComplexSolver>(comm, MatrixSymmetry::SYMMETRIC,
                                                    reorder, 0.0, true, 0);
      };
      break;
#endif
    default:
      MFEM_ABORT("Unexpected solver type!");
  }
  return std::make_unique<MeasuredComplexCoarseSolver>(
      comm,
      std::make_unique<MfemWrapperSolver<ComplexOperator>>(std::move(real_pc), false, false,
                                                           false, true),
      std::move(make_complex_pc), MeshNoise(fespace));
}

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

  std::unique_ptr<ComplexOperator> Assemble(double k, double mr, double mi) const
  {
    fem::DefaultIntegrationOrder::p_trial = order;
    auto Part = [&](double curlcurl, double mass) -> std::unique_ptr<Operator>
    {
      MaterialPropertyCoefficient df(mat_op.MaxCeedAttribute()),
          f(mat_op.MaxCeedAttribute());
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
    };
    return std::make_unique<ComplexParOperator>(Part(k, mr), Part(0.0, mi), fespace);
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

TEST_CASE("MeasuredComplexCoarseSolver chooses from the error of the real-valued solve",
          "[complexcoarse][Serial][Parallel]")
{
  // For Ar = M and Ai = c M, the real-valued solve B = (Ar + Ai)⁻¹ has the error
  // I - B A = c (1 - i) / (1 + c) I, so ‖Eᵏ w‖ / ‖w‖ = (√2 c / (1 + c))ᵏ: below 1 (the
  // real-valued approximation is kept) for c = 0.1, above 1 (the complex-valued matrix is
  // factored) for c = 10.
  ComplexCoarseChoiceReset reset;
  MPI_Comm comm = Mpi::World();
  const auto type = GENERATE(from_range(ComplexArithmeticSolvers()));
  auto serial_mesh = MakeTestMesh(3);
  TestProblem problem(comm, serial_mesh);
  const auto A_fine = problem.Assemble(0.0, 1.0, 0.1);
  const auto A_poor = problem.Assemble(0.0, 1.0, 10.0);
  auto Expected = [](double c, int k)
  { return std::pow(std::sqrt(2.0) * c / (1.0 + c), k); };

  const int n = A_fine->Height();
  ComplexVector b(n), x(n);
  b.UseDevice(true);
  x.UseDevice(true);
  b.Real().Randomize(1);
  b.Imag().Randomize(2);

  // The real-valued approximation is kept: three solves are measured, and it is applied.
  auto solver = MakeMeasuredSolver(type, comm, problem.fespace);
  solver->SetOperator(*A_fine);
  CHECK(!solver->UsesComplexMatrix());
  CHECK(!MeasuredComplexCoarseSolver::ComplexMatrixChosen());
  REQUIRE(solver->GetRealSolveError().size() == 3);
  for (int k = 0; k < 3; k++)
  {
    CHECK(solver->GetRealSolveError()[k] ==
          Catch::Approx(Expected(0.1, k + 1)).epsilon(1.0e-8));
  }
  solver->Mult(b, x);
  CHECK(Residual(comm, *A_fine, x, b) == Catch::Approx(Expected(0.1, 1)).epsilon(1.0e-8));

  // A later operator of the same solver whose real-valued approximation is poor: the
  // measurement stops at the first error above 1 and the complex-valued matrix is factored.
  solver->SetOperator(*A_poor);
  CHECK(solver->UsesComplexMatrix());
  CHECK(MeasuredComplexCoarseSolver::ComplexMatrixChosen());
  REQUIRE(solver->GetRealSolveError().size() == 1);
  CHECK(solver->GetRealSolveError()[0] == Catch::Approx(Expected(10.0, 1)).epsilon(1.0e-8));
  solver->Mult(b, x);
  CHECK(Residual(comm, *A_poor, x, b) < 1.0e-10);

  // The choice is kept for the later operators of the solver, and a later solver (as on
  // the next mesh of an adaptive refinement) factors the complex-valued matrix without
  // measuring, also where the real-valued approximation would be fine.
  solver->SetOperator(*A_fine);
  CHECK(solver->UsesComplexMatrix());
  solver->Mult(b, x);
  CHECK(Residual(comm, *A_fine, x, b) < 1.0e-10);
  auto later_solver = MakeMeasuredSolver(type, comm, problem.fespace);
  later_solver->SetOperator(*A_fine);
  CHECK(later_solver->UsesComplexMatrix());
  CHECK(later_solver->GetRealSolveError().empty());
  later_solver->Mult(b, x);
  CHECK(Residual(comm, *A_fine, x, b) < 1.0e-10);

  // Until the choice is reset (at the start of a simulation).
  MeasuredComplexCoarseSolver::ResetChoice();
  auto new_solver = MakeMeasuredSolver(type, comm, problem.fespace);
  new_solver->SetOperator(*A_fine);
  CHECK(!new_solver->UsesComplexMatrix());
  CHECK(new_solver->GetRealSolveError().size() == 3);
}

TEST_CASE("MeasuredComplexCoarseSolver does not depend on the number of processes",
          "[complexcoarse][Serial][Parallel]")
{
  // The measurement vector (MeshNoise), the measured errors and the choice on the
  // processes of the run match those of the same problem on one process, also on a
  // periodic mesh (where the elements on the two sides of a periodic boundary see different
  // coordinates for the same degree of freedom). The operators are
  // Ar = K + M, Ai = c M: the real-valued solve damps the components with a large curl and
  // amplifies the gradients by √2 c / (1 + c) per solve, so the measure decays for c = 0.2
  // and exceeds 1 at the third solve for c = 50.
  ComplexCoarseChoiceReset reset;
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
      MeasuredComplexCoarseSolver::ResetChoice();
      auto solver = MakeMeasuredSolver(type, comm, problem.fespace);
      solver->SetOperator(*problem.Assemble(1.0, 1.0, c));
      (c < 1.0 ? m.err_fine : m.err_poor) = solver->GetRealSolveError();
      (c < 1.0 ? m.complex_fine : m.complex_poor) = solver->UsesComplexMatrix();
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
