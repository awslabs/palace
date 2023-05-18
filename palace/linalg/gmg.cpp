// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "gmg.hpp"

#include "linalg/chebyshev.hpp"
#include "linalg/distrelaxation.hpp"

namespace palace
{

GeometricMultigridSolver::GeometricMultigridSolver(
    std::unique_ptr<mfem::Solver> &&coarse_solver,
    mfem::ParFiniteElementSpaceHierarchy &fespaces,
    mfem::ParFiniteElementSpaceHierarchy *aux_fespaces, int cycle_it, int smooth_it,
    int cheby_order)
  : mfem::Solver(), pc_it(cycle_it)
{
  // Configure levels of geometric coarsening. Multigrid vectors will be configured at first
  // call to Mult. The multigrid operator size is set based on the finest space dimension.
  const int n_levels = fespaces.GetNumLevels();
  MFEM_VERIFY(n_levels > 0,
              "Empty finite element space hierarchy during multigrid solver setup!");
  A_.resize(n_levels, nullptr);
  P_.resize(n_levels - 1, nullptr);
  x_.resize(n_levels, Vector());
  y_.resize(n_levels, Vector());
  r_.resize(n_levels, Vector());
  xrefs_.resize(n_levels, std::vector<Vector>());
  yrefs_.resize(n_levels, std::vector<Vector>());
  rrefs_.resize(n_levels, std::vector<Vector>());
  X_.resize(n_levels, mfem::Array<Vector *>());
  Y_.resize(n_levels, mfem::Array<Vector *>());
  R_.resize(n_levels, mfem::Array<Vector *>());

  // Configure prolongation operators.
  for (int l = 0; l < n_levels - 1; l++)
  {
    const auto *PtAP_l =
        dynamic_cast<const ParOperator *>(fespaces.GetProlongationAtLevel(l));
    MFEM_VERIFY(PtAP_l,
                "GeometricMultigridSolver requires ParOperator prolongation operators!");
    P_[l] = PtAP_l;
  }

  // Use the supplied level 0 (coarse) solver.
  B_.reserve(n_levels);
  B_.push_back(std::move(coarse_solver));

  // Configure level smoothers. Use distributive relaxation smoothing if an auxiliary
  // finite element space was provided.
  if (aux_fespaces)
  {
    for (int l = 1; l < n_levels; l++)
    {
      B_.push_back(std::make_unique<DistRelaxationSmoother>(
          fespaces.GetFESpaceAtLevel(l), aux_fespaces->GetFESpaceAtLevel(l), smooth_it, 1,
          cheby_order));
    }
  }
  else
  {
    for (int l = 1; l < n_levels; l++)
    {
      B_.push_back(std::make_unique<ChebyshevSmoother>(smooth_it, cheby_order));
    }
  }
}

void GeometricMultigridSolver::SetOperator(
    const std::vector<std::unique_ptr<ParOperator>> &ops,
    const std::vector<std::unique_ptr<ParOperator>> *aux_ops)
{
  const int n_levels = static_cast<int>(A_.size());
  MFEM_VERIFY(static_cast<std::size_t>(ops.size()) == n_levels &&
                  (!aux_ops || static_cast<std::size_t>(aux_ops->size()) == n_levels),
              "Invalid number of levels for operators in multigrid solver setup!");
  for (int l = 0; l < n_levels; l++)
  {
    A_[l] = ops[l].get();
    auto *dist_smoother = dynamic_cast<DistRelaxationSmoother *>(B_[l].get());
    if (dist_smoother)
    {
      MFEM_VERIFY(aux_ops, "Distributive relaxation smoother relies on both primary space "
                           "and auxiliary space operators for multigrid smoothing!");
      dist_smoother->SetOperator(*ops[l], *(*aux_ops)[l]);
    }
    else
    {
      B_[l]->SetOperator(*ops[l]);
    }
  }

  // Operator size is given by the fine level dimensions.
  height = A_.back()->Height();
  width = A_.back()->Width();
}

void GeometricMultigridSolver::ArrayMult(const mfem::Array<const Vector *> &X,
                                         mfem::Array<Vector *> &Y) const
{
  // Initialize.
  const int n_levels = static_cast<int>(A_.size()), n_rhs = X.Size();
  MFEM_ASSERT(!iterative_mode, "Geometric multigrid solver does not use iterative_mode!");
  MFEM_ASSERT(n_levels > 1 || pc_it == 1,
              "Single-level geometric multigrid will not work with multiple iterations!");
  if (n_rhs * height != x_.back().Size())
  {
    for (int l = 0; l < n_levels; l++)
    {
      MFEM_ASSERT(A_[l], "Missing operator for geometric multigrid level " << l << "!");
      x_[l].SetSize(n_rhs * A_[l]->Height());
      y_[l].SetSize(n_rhs * A_[l]->Height());
      r_[l].SetSize(n_rhs * A_[l]->Height());
      xrefs_[l].resize(n_rhs);
      yrefs_[l].resize(n_rhs);
      rrefs_[l].resize(n_rhs);
      X_[l].SetSize(n_rhs);
      Y_[l].SetSize(n_rhs);
      R_[l].SetSize(n_rhs);
      for (int j = 0; j < n_rhs; j++)
      {
        xrefs_[l][j].MakeRef(x_[l], j * A_[l]->Height(), A_[l]->Height());
        yrefs_[l][j].MakeRef(y_[l], j * A_[l]->Height(), A_[l]->Height());
        rrefs_[l][j].MakeRef(r_[l], j * A_[l]->Height(), A_[l]->Height());
        X_[l][j] = &xrefs_[l][j];
        Y_[l][j] = &yrefs_[l][j];
        R_[l][j] = &rrefs_[l][j];
      }
    }
  }

  // Apply V-cycle. X_ and Y_ on the finest level just point to X and Y to avoid an extra
  // copy.
  for (int j = 0; j < n_rhs; j++)
  {
    X_.back()[j] = const_cast<Vector *>(X[j]);
    Y_.back()[j] = Y[j];
  }
  for (int it = 0; it < pc_it; it++)
  {
    VCycle(n_levels - 1, (it > 0));
  }
}

void GeometricMultigridSolver::VCycle(int l, bool initial_guess) const
{
  // Pre-smooth, with zero initial guess (Y = 0 set inside). This is the coarse solve at
  // level 0. Important to note that the smoothers must respect the iterative_mode flag
  // correctly (given X, Y, compute Y <- Y + B (X - A Y)) .
  const int n_rhs = X_[l].Size();
  B_[l]->iterative_mode = initial_guess;
  B_[l]->ArrayMult(X_[l], Y_[l]);
  if (l == 0)
  {
    return;
  }

  // Compute residual.
  A_[l]->ArrayMult(Y_[l], R_[l]);
  for (int j = 0; j < n_rhs; j++)
  {
    subtract(*X_[l][j], *R_[l][j], *R_[l][j]);
  }

  // Coarse grid correction.
  P_[l - 1]->ArrayMultTranspose(R_[l], X_[l - 1]);
  if (A_[l - 1]->GetEssentialTrueDofs())
  {
    const mfem::Array<int> &dbc_tdof_list = *A_[l - 1]->GetEssentialTrueDofs();
    for (int j = 0; j < n_rhs; j++)
    {
      X_[l - 1][j]->SetSubVector(dbc_tdof_list, 0.0);
    }
  }
  VCycle(l - 1, false);

  // Prolongate and add.
  P_[l - 1]->ArrayMult(Y_[l - 1], R_[l]);
  for (int j = 0; j < n_rhs; j++)
  {
    *Y_[l][j] += *R_[l][j];
  }

  // Post-smooth, with nonzero initial guess.
  B_[l]->iterative_mode = true;
  B_[l]->ArrayMultTranspose(X_[l], Y_[l]);
}

}  // namespace palace
