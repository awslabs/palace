// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DIST_RELAXATION_SMOOTHER_HPP
#define PALACE_DIST_RELAXATION_SMOOTHER_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>

namespace palace
{

//
// Hiptmair distributive relaxation smoother applying smoothers to both the operator in the
// primary space as well as its projection into an auxiliary space.
// Reference: Hiptmair, Multigrid method for Maxwell's equations, SIAM J. Numer. Anal.
//            (1998).
//
class DistRelaxationSmoother : public mfem::Solver
{
private:
  // System matrix and its projection G^T A G (not owned).
  const mfem::Operator *A, *A_G;

  // Discrete gradient matrix.
  std::unique_ptr<mfem::Operator> G;

  // Point smoother objects for each matrix.
  mutable std::unique_ptr<mfem::Solver> B;
  std::unique_ptr<mfem::Solver> B_G;

  // Temporary vectors for smoother application.
  mutable mfem::Vector r, x_G, y_G;

  // Dirichlet boundary conditions in the auxiliary space.
  mfem::Array<int> h1_dbc_tdof_list;

  // Number of smoother iterations.
  const int pc_it;

public:
  DistRelaxationSmoother(mfem::ParFiniteElementSpace &nd_fespace,
                         mfem::ParFiniteElementSpace &h1_fespace,
                         const mfem::Array<int> &dbc_marker, int smooth_it,
                         int cheby_smooth_it, int cheby_order);

  void SetOperator(const mfem::Operator &op) override
  {
    MFEM_ABORT("SetOperator with a single operator is not implemented for "
               "DistRelaxationSmoother, use the two argument signature instead!");
  }

  void SetOperator(const mfem::Operator &op, const mfem::Operator &op_G);

  void Mult(const mfem::Vector &x, mfem::Vector &y) const override
  {
    mfem::Array<const mfem::Vector *> X(1);
    mfem::Array<mfem::Vector *> Y(1);
    X[0] = &x;
    Y[0] = &y;
    ArrayMult(X, Y);
  }

  void MultTranspose(const mfem::Vector &x, mfem::Vector &y) const override
  {
    mfem::Array<const mfem::Vector *> X(1);
    mfem::Array<mfem::Vector *> Y(1);
    X[0] = &x;
    Y[0] = &y;
    ArrayMultTranspose(X, Y);
  }

  void ArrayMult(const mfem::Array<const mfem::Vector *> &X,
                 mfem::Array<mfem::Vector *> &Y) const override
  {
    // Initialize.
    const int nrhs = X.Size();
    mfem::Array<mfem::Vector *> R(nrhs), X_G(nrhs), Y_G(nrhs);
    std::vector<mfem::Vector> rrefs(nrhs), xgrefs(nrhs), ygrefs(nrhs);
    if (nrhs * height != r.Size())
    {
      r.SetSize(nrhs * height);
      x_G.SetSize(nrhs * A_G->Height());
      y_G.SetSize(nrhs * A_G->Height());
    }
    for (int j = 0; j < nrhs; j++)
    {
      rrefs[j].MakeRef(r, j * height, height);
      xgrefs[j].MakeRef(x_G, j * A_G->Height(), A_G->Height());
      ygrefs[j].MakeRef(y_G, j * A_G->Height(), A_G->Height());
      R[j] = &rrefs[j];
      X_G[j] = &xgrefs[j];
      Y_G[j] = &ygrefs[j];
    }

    // Apply smoother.
    for (int it = 0; it < pc_it; it++)
    {
      // y = y + B (x - A y)
      B->iterative_mode = (iterative_mode || it > 0);
      B->ArrayMult(X, Y);

      // y = y + G B_G Gᵀ (x - A y)
      A->ArrayMult(Y, R);
      for (int j = 0; j < nrhs; j++)
      {
        subtract(*X[j], *R[j], *R[j]);
      }
      G->ArrayMultTranspose(R, X_G);
      for (int j = 0; j < nrhs; j++)
      {
        X_G[j]->SetSubVector(h1_dbc_tdof_list, 0.0);
      }
      B_G->ArrayMult(X_G, Y_G);
      G->ArrayAddMult(Y_G, Y, 1.0);
    }
  }

  void ArrayMultTranspose(const mfem::Array<const mfem::Vector *> &X,
                          mfem::Array<mfem::Vector *> &Y) const override
  {
    // Initialize.
    const int nrhs = X.Size();
    mfem::Array<mfem::Vector *> R(nrhs), X_G(nrhs), Y_G(nrhs);
    std::vector<mfem::Vector> rrefs(nrhs), xgrefs(nrhs), ygrefs(nrhs);
    if (nrhs * height != r.Size())
    {
      r.SetSize(nrhs * height);
      x_G.SetSize(nrhs * A_G->Height());
      y_G.SetSize(nrhs * A_G->Height());
    }
    for (int j = 0; j < nrhs; j++)
    {
      rrefs[j].MakeRef(r, j * height, height);
      xgrefs[j].MakeRef(x_G, j * A_G->Height(), A_G->Height());
      ygrefs[j].MakeRef(y_G, j * A_G->Height(), A_G->Height());
      R[j] = &rrefs[j];
      X_G[j] = &xgrefs[j];
      Y_G[j] = &ygrefs[j];
    }

    // Apply transpose.
    B->iterative_mode = true;
    for (int it = 0; it < pc_it; it++)
    {
      // y = y + G B_Gᵀ Gᵀ (x - A y)
      if (iterative_mode || it > 0)
      {
        A->ArrayMult(Y, R);
        for (int j = 0; j < nrhs; j++)
        {
          subtract(*X[j], *R[j], *R[j]);
        }
        G->ArrayMultTranspose(R, X_G);
        for (int j = 0; j < nrhs; j++)
        {
          X_G[j]->SetSubVector(h1_dbc_tdof_list, 0.0);
        }
        B_G->ArrayMultTranspose(X_G, Y_G);
        G->ArrayAddMult(Y_G, Y, 1.0);
      }
      else
      {
        G->ArrayMultTranspose(X, X_G);
        for (int j = 0; j < nrhs; j++)
        {
          X_G[j]->SetSubVector(h1_dbc_tdof_list, 0.0);
        }
        B_G->ArrayMultTranspose(X_G, Y_G);
        G->ArrayMult(Y_G, Y);
      }

      // y = y + Bᵀ (x - A y)
      B->ArrayMultTranspose(X, Y);
    }
  }
};

}  // namespace palace

#endif  // PALACE_DIST_RELAXATION_SMOOTHER_HPP
