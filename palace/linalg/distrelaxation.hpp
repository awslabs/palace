// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DIST_RELAXATION_SMOOTHER_HPP
#define PALACE_DIST_RELAXATION_SMOOTHER_HPP

#include <memory>
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
    for (int it = 0; it < pc_it; it++)
    {
      // y = y + B (x - A y)
      if (iterative_mode || it > 0)
      {
        B->iterative_mode = true;
      }
      else
      {
        B->iterative_mode = false;
      }
      B->Mult(x, y);

      // y = y + G B_G Gᵀ (x - A y)
      A->Mult(y, r);
      subtract(x, r, r);
      G->MultTranspose(r, x_G);
      x_G.SetSubVector(h1_dbc_tdof_list, 0.0);
      B_G->Mult(x_G, y_G);
      G->AddMult(y_G, y, 1.0);
    }
  }

  void MultTranspose(const mfem::Vector &x, mfem::Vector &y) const override
  {
    B->iterative_mode = true;
    for (int it = 0; it < pc_it; it++)
    {
      // y = y + G B_Gᵀ Gᵀ (x - A y)
      if (iterative_mode || it > 0)
      {
        A->Mult(y, r);
        subtract(x, r, r);
        G->MultTranspose(r, x_G);
        x_G.SetSubVector(h1_dbc_tdof_list, 0.0);
        B_G->MultTranspose(x_G, y_G);
        G->AddMult(y_G, y, 1.0);
      }
      else
      {
        G->MultTranspose(x, x_G);
        x_G.SetSubVector(h1_dbc_tdof_list, 0.0);
        B_G->MultTranspose(x_G, y_G);
        G->Mult(y_G, y);
      }

      // y = y + Bᵀ (x - A y)
      B->MultTranspose(x, y);
    }
  }
};

}  // namespace palace

#endif  // PALACE_DIST_RELAXATION_SMOOTHER_HPP
