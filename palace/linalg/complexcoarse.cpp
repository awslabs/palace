// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "complexcoarse.hpp"

#include <cmath>
#include <cstdint>
#include <utility>
#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "linalg/complexcsr.hpp"
#include "linalg/rap.hpp"
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

namespace
{

// The parallel assembled sum of a part of an operator and of its correction (ParOperator
// objects, whose parallel assembled matrices are released), or null without a correction.
std::unique_ptr<mfem::HypreParMatrix> AssembleSum(const Operator *a, const Operator *d)
{
  if (!d)
  {
    return {};
  }
  const auto *PtAPd = dynamic_cast<const ParOperator *>(d);
  MFEM_VERIFY(PtAPd, "The correction of a complex-valued coarse solve must be a "
                     "ParOperator!");
  if (!a)
  {
    return PtAPd->StealParallelAssemble();
  }
  const auto *PtAPa = dynamic_cast<const ParOperator *>(a);
  MFEM_VERIFY(PtAPa, "The operator of a complex-valued coarse solve with a correction must "
                     "be a ParOperator!");
  std::unique_ptr<mfem::HypreParMatrix> sum(
      mfem::Add(1.0, PtAPa->ParallelAssemble(), 1.0, PtAPd->ParallelAssemble()));
  PtAPa->StealParallelAssemble();
  PtAPd->StealParallelAssemble();
  return sum;
}

}  // namespace

ComplexCoarseSolver::ComplexCoarseSolver(MPI_Comm comm, SolverFactory make_real_pc,
                                         SolverFactory make_complex_pc, bool print,
                                         Vector &&w)
  : Solver<ComplexOperator>(), comm(comm), make_real_pc(std::move(make_real_pc)),
    make_complex_pc(std::move(make_complex_pc)), w(std::move(w)), print(print)
{
}

void ComplexCoarseSolver::SetInitialGuess(bool guess)
{
  Solver<ComplexOperator>::SetInitialGuess(guess);
  if (auto *pc = real_pc ? real_pc.get() : complex_pc.get())
  {
    pc->SetInitialGuess(guess);
  }
}

std::vector<double>
ComplexCoarseSolver::MeasureRealSolveError(const ComplexOperator &op,
                                           const ComplexOperator *correction) const
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
    if (correction)
    {
      correction->AddMult(u, r);
    }
    real_pc->Mult(r, z);
    u -= z;
    err.push_back(linalg::Norml2(comm, u));
    if (!(err.back() <= 1.0))
    {
      break;
    }
  }
  return err;
}

void ComplexCoarseSolver::SetOperator(const ComplexOperator &op)
{
  auto Build =
      [this](std::unique_ptr<Solver<ComplexOperator>> &pc, const SolverFactory &make_pc)
  {
    if (!pc)
    {
      pc = make_pc();
      pc->SetInitialGuess(this->initial_guess);
    }
  };
  const Matrix previous = std::exchange(matrix, Matrix::COMPLEX);
  bool *chosen = std::exchange(complex_chosen, nullptr);
  const ComplexOperator *correction = std::exchange(correction_op, nullptr);
  measured_complex = measured_complex || (chosen && *chosen);
  bool measured = false;
  if (!correction && HasZeroImag(op))
  {
    // The assembled zero imaginary part is not needed any more (the solver releases the
    // real part, as the solver in complex arithmetic releases both).
    matrix = Matrix::REAL_PART;
    complex_pc.reset();
    Build(real_pc, make_real_pc);
    if (const auto *PtAPi = dynamic_cast<const ParOperator *>(op.Imag()))
    {
      PtAPi->StealParallelAssemble();
    }
    real_pc->SetOperator(ComplexWrapperOperator(op.Real(), nullptr));
  }
  else if (w.Size() > 0 && !measured_complex)
  {
    // Factor the real-valued approximation, which is kept unless its measured error exceeds
    // 1 or is not a number.
    Build(real_pc, make_real_pc);
    real_pc->SetOperator(op);
    real_solve_error = MeasureRealSolveError(op, correction);
    measured = true;
    measured_complex = !(real_solve_error.back() <= 1.0);
    if (!measured_complex)
    {
      matrix = Matrix::REAL_APPROXIMATION;
    }
    else if (chosen)
    {
      *chosen = true;
    }
  }
  if (matrix == Matrix::COMPLEX)
  {
    // The real-valued factors are released before the complex-valued matrix is factored.
    // After a measurement, the real-valued solver has released the parallel assembled real
    // and imaginary parts of the operator, and the solver in complex arithmetic assembles
    // them again.
    real_pc.reset();
    Build(complex_pc, make_complex_pc);
    if (correction)
    {
      const auto hAr = AssembleSum(op.Real(), correction->Real()),
                 hAi = AssembleSum(op.Imag(), correction->Imag());
      complex_pc->SetOperator(
          ComplexWrapperOperator(hAr ? hAr.get() : op.Real(), hAi ? hAi.get() : op.Imag()));
    }
    else
    {
      complex_pc->SetOperator(op);
    }
  }
  if (print && matrix != previous)
  {
    Mpi::Print(comm, " Coarse solve with the {}{}\n",
               (matrix == Matrix::REAL_PART)
                   ? "real part of the complex-valued matrix in real arithmetic (its "
                     "imaginary part is zero)"
               : (matrix == Matrix::REAL_APPROXIMATION)
                   ? "real-valued approximation of the complex-valued matrix"
                   : "complex-valued matrix in complex arithmetic",
               measured ? fmt::format(" (the measured error of successive real-valued "
                                      "solves {} within 1: {:.3e})",
                                      measured_complex ? "does not stay" : "stays",
                                      fmt::join(real_solve_error, ", "))
                        : "");
  }
  this->height = op.Height();
  this->width = op.Width();
}

void ComplexCoarseSolver::Mult(const ComplexVector &x, ComplexVector &y) const
{
  (real_pc ? *real_pc : *complex_pc).Mult(x, y);
}

}  // namespace palace
