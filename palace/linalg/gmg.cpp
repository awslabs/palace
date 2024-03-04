// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "gmg.hpp"

#include <mfem.hpp>
#include "linalg/chebyshev.hpp"
#include "linalg/distrelaxation.hpp"
#include "linalg/rap.hpp"
#include "utils/timer.hpp"

namespace palace
{

template <typename OperType>
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

template <typename OperType>
void GeometricMultigridSolver<OperType>::SetOperator(const OperType &op)
{
  using ParOperType =
      typename std::conditional<std::is_same<OperType, ComplexOperator>::value,
                                ComplexParOperator, ParOperator>::type;

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

  this->height = op.Height();
  this->width = op.Width();
}

template <typename OperType>
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

template <typename OperType>
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

  // Post-smooth, with nonzero initial guess.
  B[l]->SetInitialGuess(true);
  B[l]->MultTranspose2(X[l], Y[l], R[l]);
}

template class GeometricMultigridSolver<Operator>;
template class GeometricMultigridSolver<ComplexOperator>;

}  // namespace palace
