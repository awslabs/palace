// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "gmg.hpp"

#include "linalg/chebyshev.hpp"
#include "linalg/distrelaxation.hpp"

namespace palace
{

GeometricMultigridSolver::GeometricMultigridSolver(
    std::unique_ptr<mfem::Solver> &&coarse_solver, const mfem::Array<int> &dbc_marker,
    mfem::ParFiniteElementSpaceHierarchy &fespaces,
    mfem::ParFiniteElementSpaceHierarchy *aux_fespaces, int cycle_it, int smooth_it,
    int cheby_order)
  : mfem::Solver(), fespaces_(fespaces), pc_it(cycle_it)
{
  // Read configuration file parameters used to set up the preconditioner. The default MG
  // parameters are for a V-cycle with a single pre/post smoothing iteration.
  MFEM_VERIFY(GetNumLevels() > 0,
              "Empty finite element space hierarchy during multigrid solver setup!");

  // Configure levels of geometric coarsening. Multigrid vectors will be configured at first
  // call to Mult. The multigrid operator size is set based on the finest space dimension.
  const int m = GetNumLevels();
  A_.resize(m, nullptr);
  x_.resize(m, mfem::Vector());
  y_.resize(m, mfem::Vector());
  r_.resize(m, mfem::Vector());
  X_.resize(m, mfem::Array<mfem::Vector *>());
  Y_.resize(m, mfem::Array<mfem::Vector *>());
  R_.resize(m, mfem::Array<mfem::Vector *>());

  // Use the supplied level 0 (coarse) solver.
  B_.reserve(m);
  B_.push_back(std::move(coarse_solver));

  // Configure level smoothers. Use distributive relaxation smoothing if an auxiliary
  // finite element space was provided.
  if (aux_fespaces)
  {
    int cheby_smooth_it = 1;
    for (int l = 1; l < m; l++)
    {
      B_.push_back(std::make_unique<DistRelaxationSmoother>(
          fespaces.GetFESpaceAtLevel(l), aux_fespaces->GetFESpaceAtLevel(l), dbc_marker,
          smooth_it, cheby_smooth_it, cheby_order));
    }
  }
  else
  {
    for (int l = 1; l < m; l++)
    {
      mfem::Array<int> dbc_tdof_list_l;
      fespaces.GetFESpaceAtLevel(l).GetEssentialTrueDofs(dbc_marker, dbc_tdof_list_l);
      B_.push_back(
          std::make_unique<ChebyshevSmoother>(fespaces.GetFESpaceAtLevel(l).GetComm(),
                                              dbc_tdof_list_l, smooth_it, cheby_order));
    }
  }
}

void GeometricMultigridSolver::SetOperator(
    const std::vector<std::unique_ptr<mfem::Operator>> &ops,
    const std::vector<std::unique_ptr<mfem::Operator>> *aux_ops)
{
  const int m = GetNumLevels();
  MFEM_VERIFY(ops.size() == static_cast<std::size_t>(m) &&
                  (!aux_ops || aux_ops->size() == static_cast<std::size_t>(m)),
              "Invalid number of levels for operators in multigrid solver setup!");
  for (int l = 0; l < m; l++)
  {
    A_[l] = ops[l].get();
    auto *dist_smoother = dynamic_cast<DistRelaxationSmoother *>(B_[l].get());
    if (dist_smoother)
    {
      MFEM_VERIFY(aux_ops, "Distributive relaxation smoother relies on both primary space "
                           "and auxiliary space operators for geometric multigrid!")
      dist_smoother->SetOperator(*ops[l], *(*aux_ops)[l]);
    }
    else
    {
      B_[l]->SetOperator(*ops[l]);
    }
  }
  height = A_.back()->Height();
  width = A_.back()->Width();
}

void GeometricMultigridSolver::InitVectors(int nrhs) const
{
  if (nrhs * height == x_.back().Size())
  {
    return;
  }
  DestroyVectors();
  for (int l = 0; l < GetNumLevels(); l++)
  {
    MFEM_VERIFY(A_[l], "Missing operator for geometric multigrid level " << l << "!");
    x_[l].SetSize(nrhs * A_[l]->Height());
    y_[l].SetSize(nrhs * A_[l]->Height());
    r_[l].SetSize(nrhs * A_[l]->Height());
    X_[l].SetSize(nrhs);
    Y_[l].SetSize(nrhs);
    R_[l].SetSize(nrhs);
    for (int j = 0; j < nrhs; j++)
    {
      X_[l][j] = new mfem::Vector(x_[l], j * A_[l]->Height(), A_[l]->Height());
      Y_[l][j] = new mfem::Vector(y_[l], j * A_[l]->Height(), A_[l]->Height());
      R_[l][j] = new mfem::Vector(r_[l], j * A_[l]->Height(), A_[l]->Height());
    }
  }
}

void GeometricMultigridSolver::DestroyVectors() const
{
  for (int l = 0; l < GetNumLevels(); l++)
  {
    for (int j = 0; j < X_[l].Size(); j++)
    {
      delete X_[l][j];
      delete Y_[l][j];
      delete R_[l][j];
    }
  }
}

void GeometricMultigridSolver::VCycle(int l, bool initial_guess) const
{
  // Pre-smooth, with zero initial guess (Y = 0 set inside). This is the coarse solve at
  // level 0. Important to note that the smoothers must respect the iterative_mode flag
  // correctly (given X, Y, compute Y <- Y + B (X - A Y)) .
  B_[l]->iterative_mode = initial_guess;
  B_[l]->ArrayMult(X_[l], Y_[l]);
  if (l == 0)
  {
    return;
  }

  // Compute residual and restrict.
  A_[l]->ArrayMult(Y_[l], R_[l]);
  subtract(x_[l], r_[l], r_[l]);
  GetProlongationAtLevel(l - 1).ArrayMultTranspose(R_[l], X_[l - 1]);

  // Coarse grid correction.
  VCycle(l - 1, false);

  // Prolongate and add.
  GetProlongationAtLevel(l - 1).ArrayMult(Y_[l - 1], R_[l]);
  y_[l] += r_[l];

  // Post-smooth, with nonzero initial guess.
  B_[l]->iterative_mode = true;
  B_[l]->ArrayMultTranspose(X_[l], Y_[l]);
}

}  // namespace palace
