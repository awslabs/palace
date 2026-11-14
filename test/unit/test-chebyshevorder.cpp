// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <complex>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include "fem/mesh.hpp"
#include "linalg/chebyshev.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

using json = nlohmann::json;

namespace
{

// Cube of tetrahedra with a PEC face (boundary attribute 1) and absorbing boundaries
// (attribute 2) elsewhere.
mfem::Mesh MakeCubeMesh(int n)
{
  auto smesh = mfem::Mesh::MakeCartesian3D(n, n, n, mfem::Element::TETRAHEDRON);
  for (int i = 0; i < smesh.GetNBE(); i++)
  {
    smesh.SetBdrAttribute(i, (smesh.GetBdrAttribute(i) == 1) ? 1 : 2);
  }
  smesh.SetAttributes();
  return smesh;
}

// Solve a shifted Maxwell system K + iωC - ω²M with the preconditioner of an eigenmode
// simulation of order 2 (p-multigrid with distributive relaxation smoothers), returning the
// number of Krylov iterations. The Chebyshev smoothing order is the eigenmode default
// unless smooth_order is positive, and patch is merged into the configuration.
int SolveShiftedMaxwell(const mfem::Mesh &smesh, double omega, ComplexVector &x,
                        int smooth_order = -1, const json &patch = json::object())
{
  MPI_Comm comm = Mpi::World();
  json config = {
      {"Problem", {{"Type", "Eigenmode"}, {"Verbose", 0}}},
      {"Model", {{"Mesh", "cube.mesh"}, {"L0", 1.0}}},
      {"Domains",
       {{"Materials", {{{"Attributes", {1}}, {"Permittivity", 2.0}, {"LossTan", 0.01}}}}}},
      {"Boundaries",
       {{"PEC", {{"Attributes", {1}}}},
        {"Absorbing", {{"Attributes", {2}}, {"Order", 1}}}}},
      {"Solver",
       {{"Order", 2},
        {"Eigenmode", {{"Target", 0.2}}},
        {"Linear", {{"Tol", 1.0e-10}, {"MaxIts", 400}}}}}};
  if (smooth_order > 0)
  {
    config["Solver"]["Linear"]["MGSmoothOrder"] = smooth_order;
  }
  config.merge_patch(patch);
  IoData iodata(config, false);
  iodata.model.Lc = 1.0;
  auto serial_mesh = std::make_unique<mfem::Mesh>(smesh);
  iodata.NondimensionalizeInputs(serial_mesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(comm, *serial_mesh));
  SpaceOperator space_op(iodata, mesh);
  auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  using namespace std::complex_literals;
  auto A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i, K.get(),
                                    C.get(), M.get());
  auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, 1i * omega,
                                                             -omega * omega + 0.0i, omega);
  ComplexKspSolver ksp(iodata, space_op.GetNDSpaces(), &space_op.GetH1Spaces());
  ksp.SetOperators(*A, *P);
  const int n = A->Height();
  ComplexVector b(n);
  b.UseDevice(true);
  linalg::SetRandom(comm, b);
  linalg::SetSubVector(b, space_op.GetNDDbcTDofLists().back(), 0.0);
  x.SetSize(n);
  x.UseDevice(true);
  x = 0.0;
  ksp.Mult(b, x);
  return ksp.NumTotalMultIterations();
}

// Diagonal operator with entries d, which also assembles its diagonal (for the Chebyshev
// smoother).
class DiagonalTestOperator : public Operator
{
private:
  const Vector &d;

public:
  DiagonalTestOperator(const Vector &d) : Operator(d.Size()), d(d) {}

  void Mult(const Vector &x, Vector &y) const override
  {
    y = x;
    y *= d;
  }

  void AssembleDiagonal(Vector &diag) const override { diag = d; }
};

// Check that x1 matches x0, to a tolerance far below that of the solves (the parallel
// coarse factorizations are not reproducible bit for bit).
void CheckSameSolution(MPI_Comm comm, const ComplexVector &x0, ComplexVector &x1)
{
  x1 -= x0;
  CHECK(linalg::Norml2(comm, x1) <= 1.0e-10 * linalg::Norml2(comm, x0));
}

}  // namespace

TEST_CASE("Chebyshev smoothing order raised for nearly definite operators",
          "[chebyshev][Serial][Parallel]")
{
  // With the Jacobi scaling of the smoother, a diagonal operator with positive entries is
  // definite, and the order is raised; with negative entries, it is strongly indefinite.
  MPI_Comm comm = Mpi::World();
  const int n = 200;
  for (const double lo : {1.0e-3, -0.5})
  {
    Vector d(n);
    d.UseDevice(true);
    {
      double *pd = d.HostWrite();
      const int offset = n * Mpi::Rank(comm), N = n * Mpi::Size(comm);
      for (int i = 0; i < n; i++)
      {
        pd[i] = lo + (1.0 - lo) * (offset + i + 0.5) / N;  // Never zero
      }
    }
    DiagonalTestOperator D(d);
    ChebyshevSmoother<Operator> smoother(comm, 1, 4, 1.0, 8);
    smoother.SetOperator(D);
    CHECK(smoother.GetBaseOrder() == 4);
    CHECK(smoother.GetOrder() == ((lo > 0.0) ? 8 : 4));
  }
}

TEST_CASE("Raised Chebyshev smoothing orders are kept when the GMRES probe predicts fewer "
          "iterations",
          "[chebyshev][Serial][Parallel]")
{
  // The order-2 eigenmode default raises the Chebyshev order from 4 to 8 on the nearly
  // positive semi-definite levels, and the multigrid preconditioner keeps the raised orders
  // after the GMRES probe: the solve is then the one with order 8 given explicitly.
  const auto smesh = MakeCubeMesh(4);
  ComplexVector x, x4, x8;
  const int its = SolveShiftedMaxwell(smesh, 0.5, x);
  const int its4 = SolveShiftedMaxwell(smesh, 0.5, x4, 4);
  const int its8 = SolveShiftedMaxwell(smesh, 0.5, x8, 8);
  INFO("iterations: default " << its << ", order 4 " << its4 << ", order 8 " << its8);
  CHECK(10 * its8 < 7 * its4);
  CHECK(its == its8);
  CheckSameSolution(Mpi::World(), x8, x);
}

TEST_CASE(
    "Raised Chebyshev smoothing orders are not used when the base orders converge fast",
    "[chebyshev][Serial][Parallel]")
{
  // A small two-dimensional PEC cavity at a low frequency: each solve takes a few
  // iterations, of which order 8 saves at most one, so the base order is kept.
  const auto smesh = mfem::Mesh::MakeCartesian2D(4, 4, mfem::Element::TRIANGLE);
  const json patch = {
      {"Boundaries", {{"PEC", {{"Attributes", {1, 2, 3, 4}}}}, {"Absorbing", nullptr}}}};
  ComplexVector x, x4, x8;
  const int its = SolveShiftedMaxwell(smesh, 0.3, x, -1, patch);
  const int its4 = SolveShiftedMaxwell(smesh, 0.3, x4, 4, patch);
  const int its8 = SolveShiftedMaxwell(smesh, 0.3, x8, 8, patch);
  INFO("iterations: default " << its << ", order 4 " << its4 << ", order 8 " << its8);
  CHECK(its8 < its4);
  CHECK(its == its4);
  CheckSameSolution(Mpi::World(), x4, x);
}

TEST_CASE("Raised Chebyshev smoothing orders fall back to the base orders",
          "[chebyshev][Serial][Parallel]")
{
  // With Floquet terms and a poor coarse-grid correction (a Jacobi sweep on the coarse
  // level), the order raised on the nearly positive semi-definite level slows the solve
  // down: the GMRES probe falls back to the base order, and the solve is then the one with
  // order 4 given explicitly.
  const auto smesh = MakeCubeMesh(4);
  const json patch = {
      {"Boundaries",
       {{"Periodic",
         {{"FloquetWaveVector", {2.0, 0.0, 0.0}}, {"BoundaryPairs", json::array()}}}}},
      {"Solver", {{"Linear", {{"Type", "Jacobi"}}}}}};
  ComplexVector x, x4, x8;
  const int its = SolveShiftedMaxwell(smesh, 2.0, x, -1, patch);
  const int its4 = SolveShiftedMaxwell(smesh, 2.0, x4, 4, patch);
  const int its8 = SolveShiftedMaxwell(smesh, 2.0, x8, 8, patch);
  INFO("iterations: default " << its << ", order 4 " << its4 << ", order 8 " << its8);
  CHECK(4 * its8 > 5 * its4);
  CHECK(its == its4);
  CheckSameSolution(Mpi::World(), x4, x);
}

}  // namespace palace
