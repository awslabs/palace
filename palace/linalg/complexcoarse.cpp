// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "complexcoarse.hpp"

#include <utility>
#include "linalg/complexcsr.hpp"
#include "linalg/rap.hpp"
#include "utils/communication.hpp"

namespace palace
{

ComplexCoarseSolver::ComplexCoarseSolver(MPI_Comm comm, SolverFactory make_real_pc,
                                         SolverFactory make_complex_pc, bool print)
  : Solver<ComplexOperator>(), comm(comm), make_real_pc(std::move(make_real_pc)),
    make_complex_pc(std::move(make_complex_pc)), print(print)
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

void ComplexCoarseSolver::SetOperator(const ComplexOperator &op)
{
  const bool real = HasZeroImag(op);
  auto &pc = real ? real_pc : complex_pc;
  (real ? complex_pc : real_pc).reset();
  if (!pc)
  {
    pc = real ? make_real_pc() : make_complex_pc();
    pc->SetInitialGuess(this->initial_guess);
    if (print)
    {
      Mpi::Print(comm, " Coarse solve with the {}\n",
                 real ? "real part of the complex-valued matrix in real arithmetic (its "
                        "imaginary part is zero)"
                      : "complex-valued matrix in complex arithmetic");
    }
  }
  if (real)
  {
    // The assembled zero imaginary part is not needed any more (the solver releases the
    // real part, as the solver in complex arithmetic releases both).
    if (const auto *PtAPi = dynamic_cast<const ParOperator *>(op.Imag()))
    {
      PtAPi->StealParallelAssemble();
    }
    pc->SetOperator(ComplexWrapperOperator(op.Real(), nullptr));
  }
  else
  {
    pc->SetOperator(op);
  }
  this->height = op.Height();
  this->width = op.Width();
}

void ComplexCoarseSolver::Mult(const ComplexVector &x, ComplexVector &y) const
{
  (real_pc ? *real_pc : *complex_pc).Mult(x, y);
}

}  // namespace palace
