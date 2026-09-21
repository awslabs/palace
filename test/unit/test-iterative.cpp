// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cstddef>
#include <memory>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "linalg/iterative.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace
{

using namespace Catch::Matchers;

namespace
{

// Diagonally dominant tridiagonal test operator, nonsymmetric when lower != upper. The
// apply is performed on the host, which is where the test data lives.
class TridiagonalOperator : public Operator
{
private:
  double lower, upper;

public:
  TridiagonalOperator(int n, double lower, double upper)
    : Operator(n), lower(lower), upper(upper)
  {
  }

  static double Diagonal(int i) { return 4.0 + 0.1 * i; }

  void Mult(const Vector &x, Vector &y) const override
  {
    const auto *px = x.HostRead();
    auto *py = y.HostWrite();
    for (int i = 0; i < height; i++)
    {
      py[i] = Diagonal(i) * px[i];
      if (i > 0)
      {
        py[i] += lower * px[i - 1];
      }
      if (i + 1 < height)
      {
        py[i] += upper * px[i + 1];
      }
    }
  }
};

// Jacobi preconditioner for TridiagonalOperator, needed because FGMRES requires a
// preconditioner.
class TridiagonalJacobiSolver : public Solver<Operator>
{
public:
  explicit TridiagonalJacobiSolver(int n) { height = width = n; }

  void SetOperator(const Operator &op) override {}

  void Mult(const Vector &x, Vector &y) const override
  {
    const auto *px = x.HostRead();
    auto *py = y.HostWrite();
    for (int i = 0; i < height; i++)
    {
      py[i] = px[i] / TridiagonalOperator::Diagonal(i);
    }
  }
};

// Expose the protected work vectors, to check that ReleaseWorkspace actually frees them and
// that the next solve reallocates them.
template <typename OperType>
class ProbeCgSolver : public CgSolver<OperType>
{
public:
  using CgSolver<OperType>::CgSolver;

  int WorkVectorSize() const { return this->r.Size(); }
};

template <typename OperType>
class ProbeGmresSolver : public GmresSolver<OperType>
{
public:
  using GmresSolver<OperType>::GmresSolver;

  std::size_t NumBasisVectors() const { return this->V.size(); }
};

template <typename OperType>
class ProbeFgmresSolver : public FgmresSolver<OperType>
{
public:
  using FgmresSolver<OperType>::FgmresSolver;

  std::size_t NumBasisVectors() const { return this->V.size(); }
  std::size_t NumPreconditionedBasisVectors() const { return this->Z.size(); }
};

constexpr int test_size = 24;

Vector MakeTestRhs(MPI_Comm comm)
{
  Vector b(test_size);
  b.Randomize(1 + Mpi::Rank(comm));
  return b;
}

ComplexVector MakeComplexTestRhs(MPI_Comm comm)
{
  ComplexVector b(test_size);
  b.Real().Randomize(1 + Mpi::Rank(comm));
  b.Imag().Randomize(101 + Mpi::Rank(comm));
  return b;
}

// Releasing the workspace may not change the computed solution at all, since none of the
// freed storage carries state between solves.
void CheckExactlyEqual(MPI_Comm comm, const Vector &x, const Vector &x_ref)
{
  Vector diff(x);
  diff -= x_ref;
  CHECK_THAT(linalg::Norml2(comm, diff), WithinAbs(0.0, 0.0));
}

void CheckExactlyEqual(MPI_Comm comm, const ComplexVector &x, const ComplexVector &x_ref)
{
  CheckExactlyEqual(comm, x.Real(), x_ref.Real());
  CheckExactlyEqual(comm, x.Imag(), x_ref.Imag());
}

}  // namespace

TEST_CASE("GMRES workspace release reproduces the solve", "[iterative][Serial][Parallel]")
{
  const auto comm = Mpi::World();
  TridiagonalOperator A(test_size, -1.0, -2.0);
  const auto b = MakeTestRhs(comm);

  ProbeGmresSolver<Operator> gmres(comm, 0);
  gmres.SetOperator(A);
  gmres.SetRelTol(1.0e-12);
  gmres.SetMaxIter(100);
  gmres.SetRestartDim(8);

  Vector x_ref(test_size), x(test_size);
  gmres.Mult(b, x_ref);
  const int it_ref = gmres.GetNumIterations();
  CHECK(gmres.GetConverged());
  CHECK(gmres.NumBasisVectors() == 9);

  gmres.ReleaseWorkspace();
  CHECK(gmres.NumBasisVectors() == 0);

  gmres.Mult(b, x);
  CHECK(gmres.NumBasisVectors() == 9);
  CHECK(gmres.GetNumIterations() == it_ref);
  CheckExactlyEqual(comm, x, x_ref);
}

TEST_CASE("Complex GMRES workspace release reproduces the solve",
          "[iterative][Serial][Parallel]")
{
  const auto comm = Mpi::World();
  TridiagonalOperator Ar(test_size, -1.0, -2.0), Ai(test_size, 0.25, 0.5);
  ComplexWrapperOperator A(&Ar, &Ai);
  const auto b = MakeComplexTestRhs(comm);

  ProbeGmresSolver<ComplexOperator> gmres(comm, 0);
  gmres.SetOperator(A);
  gmres.SetRelTol(1.0e-12);
  gmres.SetMaxIter(100);
  gmres.SetRestartDim(8);

  ComplexVector x_ref(test_size), x(test_size);
  gmres.Mult(b, x_ref);
  const int it_ref = gmres.GetNumIterations();
  CHECK(gmres.GetConverged());
  CHECK(gmres.NumBasisVectors() == 9);

  gmres.ReleaseWorkspace();
  CHECK(gmres.NumBasisVectors() == 0);

  gmres.Mult(b, x);
  CHECK(gmres.NumBasisVectors() == 9);
  CHECK(gmres.GetNumIterations() == it_ref);
  CheckExactlyEqual(comm, x, x_ref);
}

TEST_CASE("FGMRES workspace release reproduces the solve", "[iterative][Serial][Parallel]")
{
  const auto comm = Mpi::World();
  TridiagonalOperator A(test_size, -1.0, -2.0);
  TridiagonalJacobiSolver B(test_size);
  const auto b = MakeTestRhs(comm);

  ProbeFgmresSolver<Operator> fgmres(comm, 0);
  fgmres.SetOperator(A);
  fgmres.SetPreconditioner(B);
  fgmres.SetRelTol(1.0e-12);
  fgmres.SetMaxIter(100);
  fgmres.SetRestartDim(8);

  Vector x_ref(test_size), x(test_size);
  fgmres.Mult(b, x_ref);
  const int it_ref = fgmres.GetNumIterations();
  CHECK(fgmres.GetConverged());
  CHECK(fgmres.NumBasisVectors() == 9);
  CHECK(fgmres.NumPreconditionedBasisVectors() == 9);

  fgmres.ReleaseWorkspace();
  CHECK(fgmres.NumBasisVectors() == 0);
  CHECK(fgmres.NumPreconditionedBasisVectors() == 0);

  fgmres.Mult(b, x);
  CHECK(fgmres.NumBasisVectors() == 9);
  CHECK(fgmres.NumPreconditionedBasisVectors() == 9);
  CHECK(fgmres.GetNumIterations() == it_ref);
  CheckExactlyEqual(comm, x, x_ref);
}

TEST_CASE("CG workspace release reproduces the solve", "[iterative][Serial][Parallel]")
{
  const auto comm = Mpi::World();
  TridiagonalOperator A(test_size, -1.0, -1.0);
  const auto b = MakeTestRhs(comm);

  ProbeCgSolver<Operator> cg(comm, 0);
  cg.SetOperator(A);
  cg.SetRelTol(1.0e-12);
  cg.SetMaxIter(100);

  Vector x_ref(test_size), x(test_size);
  cg.Mult(b, x_ref);
  const int it_ref = cg.GetNumIterations();
  CHECK(cg.GetConverged());
  CHECK(cg.WorkVectorSize() == test_size);

  cg.ReleaseWorkspace();
  CHECK(cg.WorkVectorSize() == 0);

  cg.Mult(b, x);
  CHECK(cg.WorkVectorSize() == test_size);
  CHECK(cg.GetNumIterations() == it_ref);
  CheckExactlyEqual(comm, x, x_ref);
}

}  // namespace palace
