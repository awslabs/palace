// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "gmg.hpp"

#include <algorithm>
#include <vector>
#include <mfem.hpp>
#include "linalg/chebyshev.hpp"
#include "linalg/distrelaxation.hpp"
#include "linalg/rap.hpp"
#include "utils/communication.hpp"
#include "utils/timer.hpp"

namespace palace
{

template <OperatorType OperType>
GeometricMultigridSolver<OperType>::GeometricMultigridSolver(
    MPI_Comm comm, std::unique_ptr<Solver<OperType>> &&coarse_solver,
    const std::vector<const Operator *> &P, const std::vector<const Operator *> *G,
    int cycle_it, int smooth_it, int cheby_order, double cheby_sf_max, double cheby_sf_min,
    bool cheby_4th_kind)
  : Solver<OperType>(), pc_it(cycle_it), P(P.begin(), P.end()), A(P.size() + 1),
    dbc_tdof_lists(P.size()), B(P.size() + 1), X(P.size() + 1), Y(P.size() + 1),
    R(P.size() + 1), use_timer(false)
{
  // Configure levels of geometric coarsening. Multigrid vectors will be configured at first
  // call to Mult. The multigrid operator size is set based on the finest space dimension.
  const auto n_levels = P.size() + 1;
  MFEM_VERIFY(n_levels > 0,
              "Empty finite element space hierarchy during multigrid solver setup!");
  MFEM_VERIFY(!G || G->size() == n_levels,
              "Invalid input for distributive relaxation smoother auxiliary space transfer "
              "operators (mismatch in number of levels)!");

  // Use the supplied level 0 (coarse) solver.
  B[0] = std::move(coarse_solver);

  // Configure level smoothers. Use distributive relaxation smoothing if an auxiliary
  // finite element space was provided.
  for (std::size_t l = 1; l < n_levels; l++)
  {
    if (G)
    {
      const int cheby_smooth_it = 1;
      B[l] = std::make_unique<DistRelaxationSmoother<OperType>>(
          comm, *(*G)[l], smooth_it, cheby_smooth_it, cheby_order, cheby_sf_max,
          cheby_sf_min, cheby_4th_kind);
    }
    else
    {
      const int cheby_smooth_it = smooth_it;
      if (cheby_4th_kind)
      {
        B[l] = std::make_unique<ChebyshevSmoother<OperType>>(comm, cheby_smooth_it,
                                                             cheby_order, cheby_sf_max);
      }
      else
      {
        B[l] = std::make_unique<ChebyshevSmoother1stKind<OperType>>(
            comm, cheby_smooth_it, cheby_order, cheby_sf_max, cheby_sf_min);
      }
    }
  }
}

template <OperatorType OperType>
GeometricMultigridSolver<OperType>::~GeometricMultigridSolver() = default;

template <OperatorType OperType>
void GeometricMultigridSolver<OperType>::SetOperator(const OperType &op)
{
  using ParOperType = std::conditional_t<std::is_same_v<OperType, ComplexOperator>,
                                         ComplexParOperator, ParOperator>;

  const auto *mg_op = dynamic_cast<const BaseMultigridOperator<OperType> *>(&op);
  MFEM_VERIFY(mg_op, "GeometricMultigridSolver requires a MultigridOperator or "
                     "ComplexMultigridOperator argument provided to SetOperator!");

  const auto n_levels = A.size();
  MFEM_VERIFY(
      mg_op->GetNumLevels() == n_levels &&
          (!mg_op->HasAuxiliaryOperators() || mg_op->GetNumAuxiliaryLevels() == n_levels),
      "Invalid number of levels for operators in multigrid solver setup!");
  for (std::size_t l = 0; l < n_levels; l++)
  {
    A[l] = &mg_op->GetOperatorAtLevel(l);
    MFEM_VERIFY(
        A[l]->Width() == A[l]->Height() &&
            (n_levels == 1 ||
             (A[l]->Height() == ((l < n_levels - 1) ? P[l]->Width() : P[l - 1]->Height()))),
        "Invalid operator sizes for GeometricMultigridSolver!");

    const auto *PtAP_l = dynamic_cast<const ParOperType *>(&mg_op->GetOperatorAtLevel(l));
    MFEM_VERIFY(
        PtAP_l,
        "GeometricMultigridSolver requires ParOperator or ComplexParOperator operators!");
    if (l < n_levels - 1)
    {
      dbc_tdof_lists[l] = PtAP_l->GetEssentialTrueDofs();
    }

    auto *dist_smoother = dynamic_cast<DistRelaxationSmoother<OperType> *>(B[l].get());
    if (dist_smoother)
    {
      MFEM_VERIFY(mg_op->HasAuxiliaryOperators(),
                  "Distributive relaxation smoother relies on both primary space and "
                  "auxiliary space operators for multigrid smoothing!");
      dist_smoother->SetOperators(mg_op->GetOperatorAtLevel(l),
                                  mg_op->GetAuxiliaryOperatorAtLevel(l));
    }
    else
    {
      B[l]->SetOperator(mg_op->GetOperatorAtLevel(l));
    }

    X[l].SetSize(A[l]->Height());
    Y[l].SetSize(A[l]->Height());
    R[l].SetSize(A[l]->Height());
    X[l].UseDevice(true);
    Y[l].UseDevice(true);
    R[l].UseDevice(true);
  }
  SetUpPMLSubdomain(*mg_op);

  this->height = op.Height();
  this->width = op.Width();
}

template <OperatorType OperType>
void GeometricMultigridSolver<OperType>::SetUpPMLSubdomain(
    const BaseMultigridOperator<OperType> &mg_op)
{
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    if (!pml_solver)
    {
      return;
    }
    const auto &tdofs = mg_op.GetPMLTrueDofs();
    const auto *A_fine =
        dynamic_cast<const ComplexParOperator *>(&mg_op.GetFinestOperator());
    MFEM_VERIFY(A_fine, "PML subdomain correction requires a ComplexParOperator!");
    const MPI_Comm comm = A_fine->GetComm();
    HYPRE_BigInt n_pml = static_cast<HYPRE_BigInt>(tdofs.size()), n_glob = n_pml;
    Mpi::GlobalSum(1, &n_glob, comm);
    if (n_glob == 0)
    {
      pml_S.reset();
      return;
    }

    // Assemble the finest level operator (the multigrid hierarchy only uses its partially
    // assembled form, so ownership of the assembled matrices can be taken here).
    auto Assemble = [](const Operator *op) -> std::unique_ptr<mfem::HypreParMatrix>
    {
      if (!op)
      {
        return nullptr;
      }
      const auto *PtAP = dynamic_cast<const ParOperator *>(op);
      MFEM_VERIFY(PtAP, "PML subdomain correction requires ParOperator operators!");
      return PtAP->StealParallelAssemble();
    };
    auto hAr = Assemble(A_fine->Real()), hAi = Assemble(A_fine->Imag());
    const mfem::HypreParMatrix &hA = hAr ? *hAr : *hAi;

    // Selection operator S: the columns of the identity of the PML true DOFs. The PML
    // subdomain unknowns are numbered in the order of the global true DOFs and distributed
    // evenly over the processes, since the PML elements in general are not (some processes
    // may own no PML unknowns, which the sparse direct solvers do not support).
    MFEM_VERIFY(HYPRE_AssumedPartitionCheck(),
                "PML subdomain correction requires Hypre's assumed partition!");
    const int n_proc = Mpi::Size(comm), rank = Mpi::Rank(comm);
    const HYPRE_BigInt n_fine = hA.GetGlobalNumRows();
    const HYPRE_BigInt fine_start = hA.RowPart()[0];
    std::vector<HYPRE_BigInt> counts(n_proc), offsets(n_proc + 1, 0);
    MPI_Allgather(&n_pml, 1, HYPRE_MPI_BIG_INT, counts.data(), 1, HYPRE_MPI_BIG_INT, comm);
    for (int p = 0; p < n_proc; p++)
    {
      offsets[p + 1] = offsets[p] + counts[p];
    }
    auto ChunkStart = [&](int p) { return (n_glob * p) / n_proc; };
    const HYPRE_BigInt chunk_start = ChunkStart(rank), chunk_end = ChunkStart(rank + 1);

    // Send the global true DOF indices of the local PML unknowns to the processes owning
    // them in the even distribution.
    std::vector<HYPRE_BigInt> send(n_pml);
    for (HYPRE_BigInt k = 0; k < n_pml; k++)
    {
      send[k] = fine_start + tdofs[k];
    }
    auto Overlap = [](HYPRE_BigInt a0, HYPRE_BigInt a1, HYPRE_BigInt b0, HYPRE_BigInt b1)
    { return std::max<HYPRE_BigInt>(0, std::min(a1, b1) - std::max(a0, b0)); };
    std::vector<int> send_counts(n_proc), send_displs(n_proc), recv_counts(n_proc),
        recv_displs(n_proc);
    for (int p = 0; p < n_proc; p++)
    {
      send_counts[p] = static_cast<int>(
          Overlap(offsets[rank], offsets[rank + 1], ChunkStart(p), ChunkStart(p + 1)));
      send_displs[p] = static_cast<int>(
          std::clamp(ChunkStart(p), offsets[rank], offsets[rank + 1]) - offsets[rank]);
      recv_counts[p] =
          static_cast<int>(Overlap(offsets[p], offsets[p + 1], chunk_start, chunk_end));
      recv_displs[p] =
          static_cast<int>(std::clamp(offsets[p], chunk_start, chunk_end) - chunk_start);
    }
    const int m_loc = static_cast<int>(chunk_end - chunk_start);
    std::vector<HYPRE_BigInt> J(m_loc);
    MPI_Alltoallv(send.data(), send_counts.data(), send_displs.data(), HYPRE_MPI_BIG_INT,
                  J.data(), recv_counts.data(), recv_displs.data(), HYPRE_MPI_BIG_INT,
                  comm);

    // Assemble Sᵀ (one unit entry per row) and transpose.
    std::vector<int> I(m_loc + 1);
    std::vector<double> D(m_loc, 1.0);
    for (int k = 0; k <= m_loc; k++)
    {
      I[k] = k;
    }
    HYPRE_BigInt row_starts[2] = {chunk_start, chunk_end};
    HYPRE_BigInt col_starts[2] = {fine_start, hA.RowPart()[1]};
    mfem::HypreParMatrix St(comm, m_loc, n_glob, n_fine, I.data(), J.data(), D.data(),
                            row_starts, col_starts);
    pml_S.reset(St.Transpose());

    // PML subdomain operator A_PML = Sᵀ A S.
    std::unique_ptr<Operator> Ar_pml, Ai_pml;
    if (hAr)
    {
      Ar_pml.reset(mfem::RAP(hAr.get(), pml_S.get()));
    }
    if (hAi)
    {
      Ai_pml.reset(mfem::RAP(hAi.get(), pml_S.get()));
    }
    hAr.reset();
    hAi.reset();
    ComplexWrapperOperator A_pml(std::move(Ar_pml), std::move(Ai_pml));
    pml_solver->SetOperator(A_pml);
    pml_r.SetSize(pml_S->Width());
    pml_x.SetSize(pml_S->Width());
    pml_r.UseDevice(true);
    pml_x.UseDevice(true);

    Mpi::Print(" PML subdomain correction: {:d} unknowns ({:.1f}% of the finest level)\n",
               n_glob, 100.0 * static_cast<double>(n_glob) / static_cast<double>(n_fine));
  }
}

template <OperatorType OperType>
void GeometricMultigridSolver<OperType>::PMLSubdomainCorrection(int l) const
{
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    if (!pml_S || l != static_cast<int>(A.size()) - 1)
    {
      return;
    }
    BlockTimer bt(Timer::KSP_PML_SOLVE, use_timer);
    A[l]->Mult(Y[l], R[l]);
    linalg::AXPBY(1.0, X[l], -1.0, R[l]);
    pml_S->MultTranspose(R[l].Real(), pml_r.Real());
    pml_S->MultTranspose(R[l].Imag(), pml_r.Imag());
    pml_solver->Mult(pml_r, pml_x);
    pml_S->Mult(1.0, pml_x.Real(), 1.0, Y[l].Real());
    pml_S->Mult(1.0, pml_x.Imag(), 1.0, Y[l].Imag());
  }
}

template <OperatorType OperType>
void GeometricMultigridSolver<OperType>::Mult(const VecType &x, VecType &y) const
{
  // Initialize.
  const auto n_levels = A.size();
  MFEM_ASSERT(!this->initial_guess,
              "Geometric multigrid solver does not use initial guess!");
  MFEM_ASSERT(n_levels > 1 || pc_it == 1,
              "Single-level geometric multigrid will not work with multiple iterations!");

  // Apply V-cycle. The initial guess for y is zero'd at the first pre-smooth iteration.
  X.back() = x;
  for (int it = 0; it < pc_it; it++)
  {
    VCycle(n_levels - 1, (it > 0));
  }
  y = Y.back();
}

namespace
{

inline void RealMult(const Operator &op, const Vector &x, Vector &y)
{
  op.Mult(x, y);
}

inline void RealMult(const Operator &op, const ComplexVector &x, ComplexVector &y)
{
  op.Mult(x.Real(), y.Real());
  op.Mult(x.Imag(), y.Imag());
}

inline void RealMultTranspose(const Operator &op, const Vector &x, Vector &y)
{
  op.MultTranspose(x, y);
}

inline void RealMultTranspose(const Operator &op, const ComplexVector &x, ComplexVector &y)
{
  op.MultTranspose(x.Real(), y.Real());
  op.MultTranspose(x.Imag(), y.Imag());
}

}  // namespace

template <OperatorType OperType>
void GeometricMultigridSolver<OperType>::VCycle(int l, bool initial_guess) const
{
  // Pre-smooth, with zero initial guess (Y = 0 set inside). This is the coarse solve at
  // level 0. Important to note that the smoothers must respect the initial guess flag
  // correctly (given X, Y, compute Y <- Y + B (X - A Y)) .
  B[l]->SetInitialGuess(initial_guess);
  if (l == 0)
  {
    BlockTimer bt(Timer::KSP_COARSE_SOLVE, use_timer);
    B[l]->Mult(X[l], Y[l]);
    return;
  }
  B[l]->Mult2(X[l], Y[l], R[l]);

  // Compute residual.
  A[l]->Mult(Y[l], R[l]);
  linalg::AXPBY(1.0, X[l], -1.0, R[l]);

  // Coarse grid correction.
  RealMultTranspose(*P[l - 1], R[l], X[l - 1]);
  if (dbc_tdof_lists[l - 1])
  {
    linalg::SetSubVector(X[l - 1], *dbc_tdof_lists[l - 1], 0.0);
  }
  VCycle(l - 1, false);

  // Prolongate and add.
  RealMult(*P[l - 1], Y[l - 1], R[l]);
  Y[l] += R[l];

  // PML subdomain correction (finest level only), then post-smooth with nonzero initial
  // guess.
  PMLSubdomainCorrection(l);
  B[l]->SetInitialGuess(true);
  B[l]->MultTranspose2(X[l], Y[l], R[l]);
}

template class GeometricMultigridSolver<Operator>;
template class GeometricMultigridSolver<ComplexOperator>;

}  // namespace palace
