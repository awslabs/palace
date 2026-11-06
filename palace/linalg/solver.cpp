// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "solver.hpp"

#include "linalg/mumps.hpp"
#include "linalg/rap.hpp"
#include "utils/perfstats.hpp"

namespace palace
{

template <>
void MfemWrapperSolver<Operator>::SetOperator(const Operator &op)
{
  BeginSetup();
  const mfem::HypreParMatrix *Ac = nullptr;
  // Operator is always assembled as a HypreParMatrix.
  if (const auto *hA = dynamic_cast<const mfem::HypreParMatrix *>(&op))
  {
    SetInnerOperator(*hA);
    Ac = hA;
  }
  else
  {
    const auto *PtAP = dynamic_cast<const ParOperator *>(&op);
    MFEM_VERIFY(PtAP,
                "MfemWrapperSolver must be able to construct a HypreParMatrix operator!");
    if (!save_assembled)
    {
      auto hA = PtAP->StealParallelAssemble();
      SetInnerOperator(*hA);
      if (perf::Verification())
      {
        A = std::move(hA);  // Kept for the residual check
        Ac = A.get();
      }
    }
    else
    {
      Ac = &PtAP->ParallelAssemble();
      SetInnerOperator(*Ac);
    }
  }
  EndSetup(Ac);
  this->height = op.Height();
  this->width = op.Width();
}

template <>
void MfemWrapperSolver<ComplexOperator>::SetOperator(const ComplexOperator &op)
{
  BeginSetup();
  const mfem::HypreParMatrix *Ac = nullptr;
  // Assemble the real and imaginary parts, then add.
  // XX TODO: Test complex matrix assembly if coarse solve supports it.
  const mfem::HypreParMatrix *hAr = dynamic_cast<const mfem::HypreParMatrix *>(op.Real());
  const mfem::HypreParMatrix *hAi = dynamic_cast<const mfem::HypreParMatrix *>(op.Imag());
  const ParOperator *PtAPr = nullptr, *PtAPi = nullptr;
  if (op.Real() && !hAr)
  {
    PtAPr = dynamic_cast<const ParOperator *>(op.Real());
    MFEM_VERIFY(PtAPr,
                "MfemWrapperSolver must be able to construct a HypreParMatrix operator!");
    hAr = &PtAPr->ParallelAssemble();
  }
  if (op.Imag() && !hAi)
  {
    PtAPi = dynamic_cast<const ParOperator *>(op.Imag());
    MFEM_VERIFY(PtAPi,
                "MfemWrapperSolver must be able to construct a HypreParMatrix operator!");
    hAi = &PtAPi->ParallelAssemble();
  }
  if (hAr && hAi)
  {
    if (complex_matrix)
    {
      // A = [Ar, Ai]
      //     [Ai, -Ar]
      // We solve A [xr; -xi] = [br; bi]
      mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);
      mfem::Array2D<double> block_coeffs(2, 2);
      blocks(0, 0) = hAr;
      blocks(0, 1) = hAi;
      blocks(1, 0) = hAi;
      blocks(1, 1) = hAr;
      block_coeffs(0, 0) = 1.0;
      block_coeffs(0, 1) = 1.0;
      block_coeffs(1, 0) = 1.0;
      block_coeffs(1, 1) = -1.0;
      A.reset(mfem::HypreParMatrixFromBlocks(blocks, &block_coeffs));
    }
    else
    {
      // A = Ar + Ai.
      A.reset(mfem::Add(1.0, *hAr, 1.0, *hAi));
    }
    if (PtAPr)
    {
      PtAPr->StealParallelAssemble();
    }
    if (PtAPi)
    {
      PtAPi->StealParallelAssemble();
    }
    if (drop_small_entries)
    {
      DropSmallEntries();
    }
    SetInnerOperator(*A);
    if (!save_assembled && !perf::Verification())
    {
      A.reset();
    }
    Ac = A.get();
  }
  else if (hAr)
  {
    if (drop_small_entries)
    {
      A = std::make_unique<mfem::HypreParMatrix>(*hAr);
      DropSmallEntries();
      SetInnerOperator(*A);
      Ac = A.get();
    }
    else
    {
      SetInnerOperator(*hAr);
      Ac = (PtAPr && !save_assembled) ? nullptr : hAr;
    }
    if (PtAPr && !save_assembled)
    {
      PtAPr->StealParallelAssemble();
    }
  }
  else if (hAi)
  {
    if (drop_small_entries)
    {
      A = std::make_unique<mfem::HypreParMatrix>(*hAi);
      DropSmallEntries();
      SetInnerOperator(*A);
      Ac = A.get();
    }
    else
    {
      SetInnerOperator(*hAi);
      Ac = (PtAPi && !save_assembled) ? nullptr : hAi;
    }
    if (PtAPi && !save_assembled)
    {
      PtAPi->StealParallelAssemble();
    }
  }
  else
  {
    MFEM_ABORT("Empty ComplexOperator for MfemWrapperSolver!");
  }
  EndSetup(Ac);
  this->height = op.Height();
  this->width = op.Width();
}

template <>
void MfemWrapperSolver<Operator>::Mult(const Vector &x, Vector &y) const
{
  const auto idx = (setup_pending && time_setup) ? Timer::KSP_COARSE_SETUP : apply_timer;
  BlockTimer bt(idx, idx != Timer::NUM_TIMINGS);
  setup_pending = false;
  if (perf_setups >= 0)
  {
    perf::Add(perf_applies, 1);
  }
  const bool first = std::exchange(first_mult, false);
  if (first && direct)
  {
    probe.Start();
  }
  pc->Mult(x, y);
  if (first)
  {
    FinishSetup({{&x, &y}});
  }
}

template <>
void MfemWrapperSolver<ComplexOperator>::Mult(const ComplexVector &x,
                                              ComplexVector &y) const
{
  const auto idx = (setup_pending && time_setup) ? Timer::KSP_COARSE_SETUP : apply_timer;
  BlockTimer bt(idx, idx != Timer::NUM_TIMINGS);
  setup_pending = false;
  if (perf_setups >= 0)
  {
    perf::Add(perf_applies, 1);
  }
  const bool first = std::exchange(first_mult, false);
  if (first && direct)
  {
    probe.Start();
  }
  if (pc->Height() == x.Size())
  {
    mfem::Array<const Vector *> X(2);
    mfem::Array<Vector *> Y(2);
    X[0] = &x.Real();
    X[1] = &x.Imag();
    Y[0] = &y.Real();
    Y[1] = &y.Imag();
    pc->ArrayMult(X, Y);
    if (first)
    {
      FinishSetup({{X[0], Y[0]}, {X[1], Y[1]}});
    }
  }
  else
  {
    const int Nx = x.Size(), Ny = y.Size();
    Vector X(2 * Nx), Y(2 * Ny), yr, yi;
    X.UseDevice(true);
    Y.UseDevice(true);
    yr.UseDevice(true);
    yi.UseDevice(true);
    linalg::SetSubVector(X, 0, x.Real());
    linalg::SetSubVector(X, Nx, x.Imag());
    pc->Mult(X, Y);
    if (first)
    {
      FinishSetup({{&X, &Y}});
    }
    Y.ReadWrite();
    yr.MakeRef(Y, 0, Ny);
    yi.MakeRef(Y, Ny, Ny);
    // [yr; yi] is the complex conjugate of the solution
    yi *= -1.0;
    y.Real() = yr;
    y.Imag() = yi;
  }
}

template <OperatorType OperType>
void MfemWrapperSolver<OperType>::InitEagerSetup()
{
#if defined(MFEM_USE_MUMPS)
  eager_setup = (dynamic_cast<MumpsSolver *>(pc.get()) != nullptr);
#endif
}

template <OperatorType OperType>
void MfemWrapperSolver<OperType>::EnableTimer(bool setup, Timer::Index apply)
{
  time_setup = setup;
  apply_timer = apply;
  const auto role = perf::CurrentRole();
  perf_setups = perf::RegisterCounter(perf::LinearSolverPath(role, "Coarse/Setups"));
  perf_refact =
      perf::RegisterCounter(perf::LinearSolverPath(role, "Coarse/Refactorizations"));
  perf_applies = perf::RegisterCounter(perf::LinearSolverPath(role, "Coarse/Applies"));
  perf_rows = perf::RegisterCounter(perf::LinearSolverPath(role, "Coarse/Rows"),
                                    perf::Reduce::SUM, perf::Accumulate::NONE);
  perf_rows_rank = perf::RegisterCounter(perf::LinearSolverPath(role, "Coarse/RowsPerRank"),
                                         perf::Reduce::STATS, perf::Accumulate::NONE);
}

template <OperatorType OperType>
void MfemWrapperSolver<OperType>::SetInnerOperator(const mfem::HypreParMatrix &op)
{
  {
    BlockTimer bt(Timer::KSP_COARSE_SETUP, time_setup && eager_setup);
    pc->SetOperator(op);
  }
  setup_pending = true;
  if (perf_setups >= 0)
  {
    // Sparse direct solvers reuse the previous symbolic factorization when requested.
    perf::Add(perf_setups, 1);
    perf::Add(perf_refact, (direct && reorder_reuse && num_setups > 0) ? 1 : 0);
    perf::Set(perf_rows, op.Height());
    perf::Set(perf_rows_rank, op.Height());
  }
  num_setups++;
}

template <OperatorType OperType>
void MfemWrapperSolver<OperType>::BeginSetup()
{
  first_mult = true;
  A_check = nullptr;
  if (direct)
  {
    probe.Start();
  }
}

template <OperatorType OperType>
void MfemWrapperSolver<OperType>::EndSetup(const mfem::HypreParMatrix *Ac)
{
  if (direct)
  {
    probe.Stop();
  }
  // Only the main solver's coarse solve is checked: nested solvers (e.g. wave port mode
  // solves) may run on a subset of the matrix communicator.
  A_check = (perf::Verification() && role == "Main") ? Ac : nullptr;
}

template <OperatorType OperType>
void MfemWrapperSolver<OperType>::FinishSetup(
    const std::vector<std::pair<const Vector *, const Vector *>> &bx) const
{
  if (!direct)
  {
    perf::RecordHypreHierarchy(*pc, role);
    return;
  }
  probe.Stop();
  probe.Record((role == "Main") ? std::string("CoarseFactorization")
                                : "CoarseFactorization/" + role);
  // The check is collective on the matrix communicator: all ranks must agree to run it.
  int check = (A_check && A_check->Height() == bx[0].first->Size()) ? 1 : 0;
  if (A_check)
  {
    MPI_Allreduce(MPI_IN_PLACE, &check, 1, MPI_INT, MPI_MIN, A_check->GetComm());
  }
  if (check)
  {
    double res = 0.0;
    for (const auto &[b, x] : bx)
    {
      res = std::max(res, perf::RelativeResidual(*A_check, *b, *x));
    }
    const std::string base = "Verification/CoarseSolve/" + role;
    perf::SetMax(perf::RegisterCounter(base + "/MaxRelativeResidual", perf::Reduce::MAX,
                                       perf::Accumulate::MAX),
                 res);
    perf::Add(
        perf::RegisterCounter(base + "/Checks", perf::Reduce::MAX, perf::Accumulate::SUM),
        1.0);
  }
}

template <OperatorType OperType>
void MfemWrapperSolver<OperType>::DropSmallEntries()
{
  const auto nnz_before = A->NNZ();
  A->DropSmallEntries(std::pow(std::numeric_limits<double>::epsilon(), 2));
  const auto nnz_after = A->NNZ();
#if defined(MFEM_USE_MUMPS)
  if (auto *mumps = dynamic_cast<MumpsSolver *>(pc.get()))
  {
    if (reorder_reuse && (num_dropped_entries != 0) &&
        (num_dropped_entries != (nnz_before - nnz_after)))
    {
      // MUMPS errors out if there are any changes to the symmetry pattern after the first
      // factorization so we don't reuse the reordering if the number of dropped entries has
      // changed.
      mumps->SetReorderReuse(false);
    }
    else if (reorder_reuse && (num_dropped_entries == (nnz_before - nnz_after)))
    {
      // Reuse the column ordering if the number of dropped entries has not changed.
      mumps->SetReorderReuse(true);
    }
  }
#endif
  num_dropped_entries = nnz_before - nnz_after;
  Mpi::Print(" Dropping {} small entries in sparse matrix out of {} ({:.1f}%)\n",
             num_dropped_entries, nnz_before,
             static_cast<double>(num_dropped_entries) / nnz_before * 100.0);
}

template class MfemWrapperSolver<Operator>;
template class MfemWrapperSolver<ComplexOperator>;

}  // namespace palace
