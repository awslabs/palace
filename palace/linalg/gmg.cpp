// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "gmg.hpp"

#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "linalg/blockjacobi.hpp"
#include "linalg/chebyshev.hpp"
#include "linalg/distrelaxation.hpp"
#include "linalg/jacobi.hpp"
#include "linalg/rap.hpp"
#include "linalg/relaxation.hpp"
#include "linalg/vertexaggregates.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

// The smoother of a multigrid level for the given options: the Chebyshev polynomial
// smoother of the relaxation method (point Jacobi or block Jacobi), or with the discrete
// gradient G into the space of the level, the distributive relaxation smoother, which
// applies that smoother to the operator of the level and to its projection GᵀAG.
template <OperatorType OperType>
std::unique_ptr<Solver<OperType>> MakeSmoother(MPI_Comm comm,
                                               const config::SmootherData &smoother,
                                               bool print, const Operator *G = nullptr)
{
  if (G)
  {
    // The outer iterations of the distributive relaxation apply one iteration of the
    // smoother on each space.
    config::SmootherData single(smoother);
    single.smooth_it = 1;
    return std::make_unique<DistRelaxationSmoother<OperType>>(
        *G, smoother.smooth_it, MakeSmoother<OperType>(comm, single, print),
        MakeSmoother<OperType>(comm, single, print));
  }
  std::unique_ptr<Relaxation<OperType>> relaxation;
  if (smoother.block_size > 1)
  {
    relaxation = std::make_unique<BlockJacobiSmoother<OperType>>(
        comm, smoother.block_size, smoother.cheby_order, print);
  }
  else
  {
    relaxation = std::make_unique<JacobiSmoother<OperType>>(comm);
  }
  if (smoother.cheby_4th_kind)
  {
    return std::make_unique<ChebyshevSmoother<OperType>>(
        std::move(relaxation), smoother.smooth_it, smoother.cheby_order,
        smoother.cheby_sf_max);
  }
  return std::make_unique<ChebyshevSmoother1stKind<OperType>>(
      std::move(relaxation), smoother.smooth_it, smoother.cheby_order,
      smoother.cheby_sf_max, smoother.cheby_sf_min);
}

}  // namespace

config::SmootherData ResolveSmoother(const config::SmootherData &smoother,
                                     const FiniteElementSpaceHierarchy &fespaces)
{
  // A hierarchy of a single level whose coarse solve is not smoothed (a sparse direct
  // solver) has no smoother: its meshes are not examined.
  bool blocks =
      (smoother.block_size > 1 && (fespaces.GetNumLevels() > 1 || smoother.coarse));
  if (blocks)
  {
    blocks = false;
    for (std::size_t l = 0; l < fespaces.GetNumLevels() && !blocks; l++)
    {
      const auto &mesh = fespaces.GetFESpaceAtLevel(l).GetMesh();
      blocks = (mesh.Dimension() == 3 && GetVertexAggregates(mesh).any);
    }
  }
  config::SmootherData resolved(smoother);
  if (!blocks)
  {
    resolved.block_size = 1;
    resolved.coarse = false;
  }
  return resolved;
}

template <OperatorType OperType>
GeometricMultigridSolver<OperType>::GeometricMultigridSolver(
    MPI_Comm comm, std::unique_ptr<Solver<OperType>> &&coarse_solver,
    const std::vector<const Operator *> &P, const std::vector<const Operator *> *G,
    int cycle_it, const config::SmootherData &smoother, bool print)
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

  // Use the supplied level 0 (coarse) solver, with a smoother before and after it if
  // requested.
  B[0] = std::move(coarse_solver);
  if (smoother.coarse)
  {
    MFEM_VERIFY(!G || (*G)[0], "Smoothing around the coarse solver requires the discrete "
                               "gradient of the coarsest level!");
    B0 = MakeSmoother<OperType>(comm, smoother, print, G ? (*G)[0] : nullptr);
  }

  // Configure level smoothers. Use distributive relaxation smoothing if an auxiliary
  // finite element space was provided.
  for (std::size_t l = 1; l < n_levels; l++)
  {
    B[l] = MakeSmoother<OperType>(comm, smoother, print, G ? (*G)[l] : nullptr);
  }
}

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

    // The coarse solver is set first: it assembles the matrix of the coarsest level, which
    // a ParOperator then also applies for its smoother.
    if (l == 0 && B0)
    {
      B[0]->SetOperator(mg_op->GetOperatorAtLevel(0));
      Z0.SetSize(A[0]->Height());
      Z0.UseDevice(true);
    }
    auto &smoother = (l == 0 && B0) ? B0 : B[l];
    auto *dist_smoother = dynamic_cast<DistRelaxationSmoother<OperType> *>(smoother.get());
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
      smoother->SetOperator(mg_op->GetOperatorAtLevel(l));
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
  if (l == 0 && B0)
  {
    // Pre-smooth, solve for the residual, and post-smooth.
    B0->SetInitialGuess(initial_guess);
    B0->Mult2(X[0], Y[0], R[0]);
    A[0]->Mult(Y[0], R[0]);
    linalg::AXPBY(1.0, X[0], -1.0, R[0]);
    {
      BlockTimer bt(Timer::KSP_COARSE_SOLVE, use_timer);
      B[0]->SetInitialGuess(false);
      B[0]->Mult(R[0], Z0);
    }
    Y[0] += Z0;
    B0->SetInitialGuess(true);
    B0->MultTranspose2(X[0], Y[0], R[0]);
    return;
  }
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
