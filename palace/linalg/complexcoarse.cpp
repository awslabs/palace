// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "complexcoarse.hpp"

#include <cmath>
#include <cstdint>
#include <utility>
#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "utils/communication.hpp"

namespace palace
{

Vector MeshNoise(const FiniteElementSpace &fespace)
{
  auto Mix = [](std::uint64_t h)
  {
    // SplitMix64.
    h += 0x9e3779b97f4a7c15ull;
    h = (h ^ (h >> 30)) * 0xbf58476d1ce4e5b9ull;
    h = (h ^ (h >> 27)) * 0x94d049bb133111ebull;
    return h ^ (h >> 31);
  };
  const auto &unbroken = fespace.GetUnbroken();
  const bool vector = (unbroken.GetFEColl().GetRangeType(unbroken.Dimension()) ==
                       mfem::FiniteElement::VECTOR);
  mfem::VectorFunctionCoefficient noise(
      vector ? unbroken.SpaceDimension() : unbroken.GetVDim(),
      [&Mix](const mfem::Vector &x, mfem::Vector &f)
      {
        std::uint64_t h = 0;
        for (int d = 0; d < x.Size(); d++)
        {
          h = Mix(h ^ static_cast<std::uint64_t>(std::llround(std::ldexp(x(d), 26))));
        }
        for (int d = 0; d < f.Size(); d++)
        {
          h = Mix(h);
          f(d) = std::ldexp(static_cast<double>(h >> 11), -52) - 1.0;
        }
      });
  // Averaged over the elements of each degree of freedom: the elements on the two sides of
  // a periodic boundary see different coordinates, and the plain projection would take the
  // value of whichever element comes last on the process which owns the degree of freedom.
  mfem::ParGridFunction gf(const_cast<mfem::ParFiniteElementSpace *>(&unbroken.Get()));
  gf.ProjectDiscCoefficient(noise, mfem::GridFunction::ARITHMETIC);
  Vector wu(unbroken.GetTrueVSize());
  gf.GetTrueDofs(wu);
  Vector w(fespace.GetTrueVSize());
  const auto *WU = wu.HostRead();
  auto *W = w.HostWrite();
  for (int t = 0; t < wu.Size(); t++)
  {
    for (int v = 0; v < fespace.GetNumVersions(t); v++)
    {
      W[fespace.GetVersionTrueDof(t, v)] = (WU[t] < 0.0) ? -1.0 : 1.0;
    }
  }
  w.UseDevice(true);
  return w;
}

MeasuredComplexCoarseSolver::MeasuredComplexCoarseSolver(
    MPI_Comm comm, std::unique_ptr<Solver<ComplexOperator>> &&real_pc,
    std::unique_ptr<Solver<ComplexOperator>> &&complex_pc, Vector &&w)
  : Solver<ComplexOperator>(), comm(comm), real_pc(std::move(real_pc)),
    complex_pc(std::move(complex_pc)), w(std::move(w))
{
}

void MeasuredComplexCoarseSolver::SetInitialGuess(bool guess)
{
  Solver<ComplexOperator>::SetInitialGuess(guess);
  if (real_pc)
  {
    real_pc->SetInitialGuess(guess);
  }
  complex_pc->SetInitialGuess(guess);
}

std::vector<double>
MeasuredComplexCoarseSolver::MeasureRealSolveError(const ComplexOperator &op) const
{
  const int n = op.Height();
  MFEM_VERIFY(w.Size() == n, "Invalid size of the measurement vector of the complex-valued "
                             "coarse solve!");
  ComplexVector u(n), r(n), z(n);
  u.UseDevice(true);
  r.UseDevice(true);
  z.UseDevice(true);
  u.Real() = w;
  u.Imag() = 0.0;
  linalg::Normalize(comm, u);
  std::vector<double> err;
  for (int k = 0; k < 3; k++)
  {
    // u ← (I - B A) u.
    op.Mult(u, r);
    real_pc->Mult(r, z);
    u -= z;
    err.push_back(linalg::Norml2(comm, u));
    if (!(err.back() > 0.0 && err.back() <= 1.0))
    {
      break;
    }
  }
  return err;
}

void MeasuredComplexCoarseSolver::SetOperator(const ComplexOperator &op)
{
  if (real_pc)
  {
    real_pc->SetOperator(op);
    if (op.Real() && op.Imag())
    {
      real_solve_error = MeasureRealSolveError(op);
      if (!(real_solve_error.back() <= 1.0))
      {
        // Also for a value which is not a number (a failed real-valued factorization). The
        // real-valued factors are released before the complex-valued matrix is factored.
        real_pc.reset();
      }
      if (!printed || !real_pc)
      {
        Mpi::Print(comm,
                   " Coarse solve with the {} (error of successive real-valued solves: "
                   "{:.3e})\n",
                   real_pc ? "real-valued approximation" : "complex-valued matrix",
                   fmt::join(real_solve_error, ", "));
        printed = true;
      }
    }
  }
  if (!real_pc)
  {
    // The real-valued solver has released the parallel assembled real and imaginary parts
    // of an operator which it was given, and the complex-valued solver assembles them
    // again.
    complex_pc->SetOperator(op);
  }
  this->height = op.Height();
  this->width = op.Width();
}

void MeasuredComplexCoarseSolver::Mult(const ComplexVector &x, ComplexVector &y) const
{
  (real_pc ? *real_pc : *complex_pc).Mult(x, y);
}

}  // namespace palace
